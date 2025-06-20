/*------------------------------------------------------------------------------
* cssr.c : Compact SSR message decode functions
*
*          Copyright (C) 2015- by Mitsubishi Electric Corp., All rights reserved.
*
* references :
*     see rtcm3.c
*
* version : $Revision:$ $Date:$
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "cssr.h"
#define CLASGRID_GLOBAL_DEFINE
#include "clasgrid.h"

#define L6FRMPREAMB 0x1ACFFC1Du /* L6 message frame preamble */
#define BLEN_MSG 218

cssr_t _cssr = {0,};
#define L6_CH_NUM 2

enum {
    ref_mask = 0,
    ref_orbit,
    ref_clock,
    ref_cbias,
    ref_pbias,
    ref_bias,
    ref_ura,
    ref_stec,
    ref_grid,
    ref_service,
    ref_combined,
    ref_atmospheric
};

/* constants -----------------------------------------------------------------*/

/* ssr update intervals ------------------------------------------------------*/
static const double ssrudint[16]={
    1,2,5,10,15,30,60,120,240,300,600,900,1800,3600,7200,10800
};

static int l6delivery[L6_CH_NUM] = {-1, -1};
static int l6facility[L6_CH_NUM] = {-1, -1};

#define MAX_NGRID   4           /* number of grids for interpolation */
#define MAX_DIST    100.0       /* max distance to grid (km) */
#define MAX_AGE     300.0       /* max age of difference (s) */

#define SQR(x)      ((x)*(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))

static double decode_sval(unsigned char *buff, int i, int n, double lsb)
{
    int slim=-((1<<(n-1))-1)-1,v;
    v = getbits(buff, i, n);
    return (v==slim) ? INVALID_VALUE:(double)v*lsb;
}

static int sys2gnss(int sys, int *prn_min)
{
    int id = CSSR_SYS_NONE;

    if (prn_min) {
        *prn_min = 1;
    }

    switch (sys) {
        case SYS_GPS: id = CSSR_SYS_GPS; break;
        case SYS_GLO: id = CSSR_SYS_GLO; break;
        case SYS_GAL: id = CSSR_SYS_GAL; break;
        case SYS_CMP: id = CSSR_SYS_BDS; break;
        case SYS_SBS:
            id = CSSR_SYS_SBS;
            if (prn_min) {
                *prn_min = 120;
            }
            break;
        case SYS_QZS:
            id = CSSR_SYS_QZS;
            if (prn_min) {
                *prn_min = 193;
            }
            break;
    }

    return id;
}

/* convert GNSS ID of cssr to system id of rtklib */
static int gnss2sys(int id, int *prn_min)
{
    int sys = SYS_NONE;

    if (prn_min) {
        *prn_min = 1;
    }

    switch (id) {
        case CSSR_SYS_GPS: sys = SYS_GPS; break;
        case CSSR_SYS_GLO: sys = SYS_GLO; break;
        case CSSR_SYS_GAL: sys = SYS_GAL; break;
        case CSSR_SYS_BDS: sys = SYS_CMP; break;
        case CSSR_SYS_SBS:
            sys = SYS_SBS;
            if (prn_min) {
                *prn_min = 120;
            }
            break;
        case CSSR_SYS_QZS:
            sys = SYS_QZS;
            if (prn_min) {
                *prn_min = 193;
            }
            break;
    }

    return sys;
}

/*
 * count number of satellite in satellite mask
 */
static int svmask2nsat(uint64_t svmask)
{
    int j,nsat=0;

    for (j=0;j<CSSR_MAX_SV_GNSS;j++) {
        if ((svmask>>(CSSR_MAX_SV_GNSS-1-j))&1) {
            nsat++;
        }
    }

    return nsat;
}

/*
 * count number of signals in signal mask
 */
static int sigmask2nsig(uint16_t sigmask)
{
    int j,nsig=0;

    for (j=0;j<CSSR_MAX_SIG;j++) {
        if ((sigmask>>j)&1) {
            nsig++;
        }
    }
    return nsig;
}


static int svmask2nsatlist(uint64_t svmask, int id, int *sat)
{
    int j,nsat=0,sys,prn_min;

    sys = gnss2sys(id, &prn_min);
    for (j=0;j<CSSR_MAX_SV_GNSS;j++) {
        if ((svmask>>(CSSR_MAX_SV_GNSS-1-j))&1) {
            sat[nsat++] = satno(sys, prn_min+j);
        }
    }
    return nsat;
}

/* convert from svmask to satellite list */
static int svmask2sat(uint64_t *svmask,int *sat)
{
    int j,id,nsat=0,sys,prn_min;

    for (id=0;id<CSSR_MAX_GNSS;id++) {
        sys = gnss2sys(id, &prn_min);
        for (j=0;j<CSSR_MAX_SV_GNSS;j++) {
            if ((svmask[id]>>(CSSR_MAX_SV_GNSS-1-j))&1) {
                if (sat)
                    sat[nsat] = satno(sys, prn_min+j);
                nsat++;
            }
        }
    }
    return nsat;
}

/* decode stec quality indicator */
static float decode_cssr_quality_stec(int a, int b)
{
    float quality;

    if ((a == 0 && b == 0) || (a == 7 && b == 7)) {
        quality = 9999 * 1000;
    } else {
        quality = (1.0+b*0.25)*pow(3.0,a)-1.0;
    }

    return quality;
}

/* decode tropo quality indicator */
static float decode_cssr_quality_trop(int a, int b)
{
    float quality;

    if ((a == 0 && b == 0) || (a == 7 && b == 7)) {
        quality = 9999;
    } else {
        quality = (1.0+b*0.25)*pow(3.0,a)-1.0;
    }

    return quality;
}

static void check_week_ref(rtcm_t *rtcm, int tow, int i)
{
    if (rtcm->obs_ref[i].time != 0 || rtcm->obs_ref[i].sec != 0.0) {
        gtime_t time = gpst2time(rtcm->week_ref[i], tow);
        if (timediff(time, rtcm->obs_ref[i]) > (86400*7/2)) {
            char temp1[64], temp2[64];
            time2str(rtcm->obs_ref[i], temp1, 0);
            time2str(time,             temp2, 0);
            trace(2, "check_week_ref(): CSSR message time is big, subtype=%2d, time=%s, obstime=%s\n",
                i + 1, temp2, temp1);
            trace(2, "check_week_ref(): adjust reference week, subtype=%2d, week=%d\n",
                i + 1, --rtcm->week_ref[i]);
        }
        rtcm->obs_ref[i].time = rtcm->obs_ref[i].sec = 0.0;
    }
    if (rtcm->tow0 != -1) {
        if (rtcm->tow_ref[i] != -1 && ((tow - rtcm->tow_ref[i]) < (-86400*7/2))) {
            ++rtcm->week_ref[i];
        }
        rtcm->tow_ref[i] = tow;
    }
}

typedef struct _cssr_bank {
    /* iono */
    stecd_t stecdata[RTCM_SSR_MAX_GP][MAXSAT];
    stec_t stec[RTCM_SSR_MAX_GP];
    /* trop */
    zwdd_t zwddata[RTCM_SSR_MAX_GP];
    zwd_t zwd[RTCM_SSR_MAX_GP];
    /* bias */
    double cbias[MAXSAT][MAXCODE];
    double pbias[MAXSAT][MAXCODE];
    int smode[MAXSAT][MAXCODE];
    /* orbit */
    double deph0[MAXSAT];
    double deph1[MAXSAT];
    double deph2[MAXSAT];
    int iode[MAXSAT];
    /* clock */
    double c0[MAXSAT];
    /* other */
    gtime_t time[MAXSAT][6];
    double udi[MAXSAT][6];
    int iod[MAXSAT][6];
    int prn[MAXSAT][6];
    int flag[MAXSAT];
    gtime_t update_time;
    gtime_t orbit_time;
    gtime_t clock_time;
    gtime_t bias_time;
    gtime_t trop_time;
    int facility;
    int gridnum;
    int network;
    int use;
} cssr_bank;

typedef struct _cssr_orbit_bank {
    int use;
    gtime_t time;
    int network;
    int prn[MAXSAT];
    double udi[MAXSAT];
    int iod[MAXSAT];
    int iode[MAXSAT];
    double deph0[MAXSAT];
    double deph1[MAXSAT];
    double deph2[MAXSAT];
} cssr_orbit_bank;

typedef struct _cssr_bias_bank {
    int use;
    gtime_t time;
    int network;
    int bflag;
    int prn[MAXSAT];
    double udi[MAXSAT];
    int iod[MAXSAT];
    int smode[MAXSAT][MAXCODE];
    int sflag[MAXSAT][MAXCODE];
    double cbias[MAXSAT][MAXCODE];
    double pbias[MAXSAT][MAXCODE];
} cssr_bias_bank;

typedef struct _cssr_clock_bank {
    int use;
    gtime_t time;
    int network;
    int prn[MAXSAT];
    double udi[MAXSAT];
    int iod[MAXSAT];
    double c0[MAXSAT];
} cssr_clock_bank;

typedef struct _cssr_trop_bank {
    int use;
    gtime_t time;
    int gridnum[CSSR_MAX_NETWORK];
    double gridpos[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][3];
    double total[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    double wet[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    int satnum[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    int prn[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    double iono[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
} cssr_trop_bank;

typedef struct _cssr_latest_trop {
    int gridnum[CSSR_MAX_NETWORK];
    double gridpos[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][3];
    gtime_t troptime[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    double total[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    double wet[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    gtime_t stectime[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    int prn[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    double stec0[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    double stec[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
} cssr_latest_trop;

#define BIASBANKNUM  128
#define TROPBANKNUM  128
#define ORBBANKNUM   128
#define CLKBANKNUM   128

typedef struct cssr_bank_control {
    cssr_orbit_bank OrbitBank[ORBBANKNUM];
    cssr_clock_bank ClockBank[CLKBANKNUM];
    cssr_bias_bank  BiasBank[BIASBANKNUM];
    cssr_trop_bank  TropBank[TROPBANKNUM];
    int fastfix[CSSR_MAX_NETWORK];
    cssr_latest_trop LatestTrop;
    int Facility;
    int separation;
    int NextOrbit;
    int NextClock;
    int NextBias;
    int NextTrop;
    int use;
} cssr_bank_control;

static cssr_bank_control cssrObject[L6_CH_NUM]; /* cssr_bank_control object */
static int chidx = 0; /* current reading cssr object idx */
grid_t BackupGrid[L6_CH_NUM];
cssr_bank CurrentCSSR[L6_CH_NUM];
cssr_bank BackupCSSR[L6_CH_NUM];

static void check_cssr_changed_facility(int facility)
{
    if (cssrObject[chidx].Facility != facility) {
        trace(4, "CCSR bank clear: facility changed, %d ---> %d\n", cssrObject[chidx].Facility + 1, facility + 1);
        memset(&cssrObject[chidx].LatestTrop, 0x00, sizeof(cssr_latest_trop));
        memset(cssrObject[chidx].OrbitBank, 0x00, sizeof(cssrObject[chidx].OrbitBank));
        memset(cssrObject[chidx].ClockBank, 0x00, sizeof(cssrObject[chidx].ClockBank));
        memset(cssrObject[chidx].BiasBank,  0x00, sizeof(cssrObject[chidx].BiasBank));
        memset(cssrObject[chidx].TropBank,  0x00, sizeof(cssrObject[chidx].TropBank));
        cssrObject[chidx].Facility = facility;
        cssrObject[chidx].separation = 0;
        cssrObject[chidx].NextOrbit = 0;
        cssrObject[chidx].NextClock = 0;
        cssrObject[chidx].NextBias  = 0;
        cssrObject[chidx].NextTrop  = 0;
    }
}

static void clear_cssr_bank(int tow, int iod)
{
    trace(1, "change cssr iod: tow=%d, iod=%d\n", time, iod);
    memset(&cssrObject[chidx].LatestTrop, 0x00, sizeof(cssr_latest_trop));
    memset(cssrObject[chidx].OrbitBank, 0x00, sizeof(cssrObject[chidx].OrbitBank));
    memset(cssrObject[chidx].ClockBank, 0x00, sizeof(cssrObject[chidx].ClockBank));
    memset(cssrObject[chidx].BiasBank,  0x00, sizeof(cssrObject[chidx].BiasBank));
    memset(cssrObject[chidx].TropBank,  0x00, sizeof(cssrObject[chidx].TropBank));
    cssrObject[chidx].separation = 0;
    cssrObject[chidx].NextOrbit = 0;
    cssrObject[chidx].NextClock = 0;
    cssrObject[chidx].NextBias  = 0;
    cssrObject[chidx].NextTrop  = 0;
}

static cssr_orbit_bank *get_same_orbit_correction(gtime_t time, int network)
{
    int i;

    for (i = 0; i < ORBBANKNUM; ++i) {
        if (cssrObject[chidx].OrbitBank[i].use == TRUE && cssrObject[chidx].OrbitBank[i].network == network &&
            timediff(cssrObject[chidx].OrbitBank[i].time, time) == 0.0) {
            return &cssrObject[chidx].OrbitBank[i];
        }
    }
    return NULL;
}

static cssr_orbit_bank *get_close_orbit_correction(gtime_t time, int network, double age)
{
    int search = ((cssrObject[chidx].separation & (1 << (network - 1))) ? network: 0);
    int pos = -1;
    int i;

    for (i = 0; i < ORBBANKNUM; ++i) {
        if (cssrObject[chidx].OrbitBank[i].use == TRUE) {
            trace(5, "get_close_orbit_correction(): orbit=%.1f, network=%d, diff=%.1f\n", time2gpst(cssrObject[chidx].OrbitBank[i].time, NULL),
                cssrObject[chidx].OrbitBank[i].network, timediff(cssrObject[chidx].OrbitBank[i].time, time));
        }
        if (cssrObject[chidx].OrbitBank[i].use == TRUE && cssrObject[chidx].OrbitBank[i].network == search && fabs(timediff(cssrObject[chidx].OrbitBank[i].time, time)) <= age) {
            if (pos != -1 && fabs(timediff(cssrObject[chidx].OrbitBank[pos].time, time)) < fabs(timediff(cssrObject[chidx].OrbitBank[i].time, time))) {
                continue;
            }
            if (timediff(cssrObject[chidx].OrbitBank[i].time, time) > 0.0) {
                continue;
            }
            pos = i;
        }
    }
    if (pos != -1) {
        trace(4, "get_close_orbit_correction(): orbit=%.1f, network=%d, diff=%.1f\n", time2gpst(cssrObject[chidx].OrbitBank[pos].time, NULL),
            cssrObject[chidx].OrbitBank[pos].network, timediff(cssrObject[chidx].OrbitBank[pos].time, time));
        return &cssrObject[chidx].OrbitBank[pos];
    }
    return NULL;
}

void set_cssr_bank_orbit(gtime_t time, nav_t *nav, int network)
{
    cssr_orbit_bank *orbit;
    int i;

    if ((orbit = get_same_orbit_correction(time, network)) == NULL) {
        cssrObject[chidx].OrbitBank[cssrObject[chidx].NextOrbit].network = network;
        cssrObject[chidx].OrbitBank[cssrObject[chidx].NextOrbit].time = time;
        cssrObject[chidx].OrbitBank[cssrObject[chidx].NextOrbit].use = TRUE;
        orbit = &cssrObject[chidx].OrbitBank[cssrObject[chidx].NextOrbit];
        if (++cssrObject[chidx].NextOrbit >= ORBBANKNUM) {
            cssrObject[chidx].NextOrbit = 0;
        }
        memset(orbit->prn, 0x00, sizeof(orbit->prn));
    }

    for (i = 0; i < MAXSAT; ++i) {
        if (orbit->prn[i] == 0) {
            continue;
        }
        orbit->prn[i] = 0;
        orbit->udi[i] = 0.0;
        orbit->iod[i] = 0;
        orbit->iode[i] = 0;
        orbit->deph0[i] = 0.0;
        orbit->deph1[i] = 0.0;
        orbit->deph2[i] = 0.0;
    }

    for (i = 0; i < MAXSAT; ++i) {
        if (nav->ssr[i].update_oc == 1) {
            if (nav->ssr[i].deph[0] != INVALID_VALUE) {
                trace(3, "set_cssr_bank_orbit(): tow=%.1f, network=%d, prn=%2d, iode=%d, deph=%.4f\n", time2gpst(time, NULL),
                    network, i + 1, nav->ssr[i].iode, nav->ssr[i].deph[0]);
                orbit->prn[i] = i + 1;
                orbit->udi[i] = nav->ssr[i].udi[0];
                orbit->iod[i] = nav->ssr[i].iod[0];
                orbit->iode[i] = nav->ssr[i].iode;
                orbit->deph0[i] = nav->ssr[i].deph[0];
                orbit->deph1[i] = nav->ssr[i].deph[1];
                orbit->deph2[i] = nav->ssr[i].deph[2];
            } else {
                trace(3, "set_cssr_bank_orbit(): tow=%.1f, network=%d, prn=%2d, iode=%d, deph=#N/A\n", time2gpst(time, NULL),
                    network, i + 1, nav->ssr[i].iode);
                orbit->prn[i] = 0;
            }
        }
    }
    trace(4, "set_cssr_bank_orbit(): tow=%d, network=%d, next=%d\n",
        (int)time2gpst(time, NULL), network, cssrObject[chidx].NextOrbit);
}

static cssr_bias_bank *get_same_bias_correction(gtime_t time, int network)
{
    int i;

    for (i = 0; i < BIASBANKNUM; ++i) {
        if (cssrObject[chidx].BiasBank[i].use == TRUE && cssrObject[chidx].BiasBank[i].network == network &&
            timediff(cssrObject[chidx].BiasBank[i].time, time) == 0.0) {
            return &cssrObject[chidx].BiasBank[i];
        }
    }
    return NULL;
}

static cssr_bias_bank *get_close_cbias_correction(gtime_t time, int network, double age)
{
    int pos = -1;
    int i;

    trace(3, "get_close_cbias_correction(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);

    if (network >= 0 && network < CSSR_MAX_NETWORK) {
        for (i = 0; i < BIASBANKNUM; ++i) {
            if (cssrObject[chidx].BiasBank[i].use == TRUE) {
                trace(5, "get_close_cbias_correction(): bias=%.1f, network=%d, flag=0x%02x\n", time2gpst(cssrObject[chidx].BiasBank[i].time, NULL),
                    cssrObject[chidx].BiasBank[i].network, cssrObject[chidx].BiasBank[i].bflag);
            }
            if (cssrObject[chidx].BiasBank[i].use == TRUE && cssrObject[chidx].BiasBank[i].network == network && timediff(cssrObject[chidx].BiasBank[i].time, time) <= age) {
                if (pos != -1 && timediff(cssrObject[chidx].BiasBank[pos].time, cssrObject[chidx].BiasBank[i].time) > 0.0) {
                    continue;
                }
                if ((cssrObject[chidx].BiasBank[i].bflag & 0x01) == 0x01) {
                    pos = i;
                }
            }
        }
        if (pos != -1) {
            trace(4, "get_close_cbias_correction(): bias=%.1f, diff=%.1f, network=%d, flag=0x%02x\n",
                time2gpst(cssrObject[chidx].BiasBank[pos].time, NULL), timediff(cssrObject[chidx].BiasBank[pos].time, time),
                cssrObject[chidx].BiasBank[pos].network, cssrObject[chidx].BiasBank[pos].bflag);
            return &cssrObject[chidx].BiasBank[pos];
        }
    }
    return NULL;
}

static cssr_bias_bank *get_close_pbias_correction(gtime_t time, int network, double age)
{
    int pos = -1;
    int i;

    trace(3, "get_close_pbias_correction(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);

    if (network >= 0 && network < CSSR_MAX_NETWORK) {
        for (i = 0; i < BIASBANKNUM; ++i) {
            if (cssrObject[chidx].BiasBank[i].use == TRUE) {
                trace(5, "get_close_pbias_correction(): bias=%.1f, network=%d, flag=0x%02x\n", time2gpst(cssrObject[chidx].BiasBank[i].time, NULL),
                    cssrObject[chidx].BiasBank[i].network, cssrObject[chidx].BiasBank[i].bflag);
            }
            if (cssrObject[chidx].BiasBank[i].use == TRUE && cssrObject[chidx].BiasBank[i].network == network && timediff(cssrObject[chidx].BiasBank[i].time, time) <= age) {
                if (pos != -1 && timediff(cssrObject[chidx].BiasBank[pos].time, cssrObject[chidx].BiasBank[i].time) > 0.0) {
                    continue;
                }
                if ((cssrObject[chidx].BiasBank[i].bflag & 0x02) == 0x02) {
                    pos = i;
                }
            }
        }
        if (pos != -1) {
            trace(4, "get_close_pbias_correction(): bias=%.1f, diff=%.1f, network=%d, flag=0x%02x\n",
                time2gpst(cssrObject[chidx].BiasBank[pos].time, NULL), timediff(cssrObject[chidx].BiasBank[pos].time, time),
                cssrObject[chidx].BiasBank[pos].network, cssrObject[chidx].BiasBank[pos].bflag);
            return &cssrObject[chidx].BiasBank[pos];
        }
    }
    return NULL;
}

static int get_bias_save_point(cssr_bias_bank *bias, int prn, int mode)
{
    int i;

    for (i = 0; i < MAXCODE; ++i) {
        if (bias->smode[prn-1][i] == mode) {
            return i;
        }
    }
    for (i = 0; i < MAXCODE; ++i) {
        if (bias->smode[prn-1][i] == 0) {
            return i;
        }
    }
    return -1;
}

static cssr_bias_bank *add_base_value_to_cbias(cssr_bias_bank *srcbias, int network)
{
    gtime_t basetime = timeadd(srcbias->time, fmod(time2gpst(srcbias->time, NULL), 30.0) * -1.0);
    static cssr_bias_bank cbias;
    cssr_bias_bank *temp;
    int i, j, pos;
    
    if ((temp = get_same_bias_correction(basetime, network)) && (temp->bflag & 0x01)) {
        trace(4, "add_base_value_to_cbias(): bias_tow=%d, network=%d, base_tow=%d, network=%d, flag=0x%02x\n", (int)time2gpst(srcbias->time, NULL),
            srcbias->network, (int)time2gpst(temp->time, NULL), temp->network, temp->bflag);
        memcpy(&cbias, srcbias, sizeof(cssr_bias_bank));
        
        for (i = 0; i < MAXSAT; ++i) {
            if (cbias.prn[i] != 0 && temp->prn[i] != 0) {
                for (j = 0; j < MAXCODE; ++j) {
                    if (temp->smode[i][j] != 0 && (temp->sflag[i][j] & 0x01)) {
                        trace(5, "add_base_value_to_cbias(): j=%2d, prn=%2d, tow=%d, mode=%2d, flag=0x%02x, value=%.3f\n",
                            j, temp->prn[i], (int)time2gpst(temp->time, NULL), temp->smode[i][j],
                            temp->sflag[i][j], temp->cbias[i][j]);
                        if ((pos = get_bias_save_point(&cbias, cbias.prn[i], temp->smode[i][j])) != -1) {
                            if (cbias.cbias[i][pos] != INVALID_VALUE && temp->cbias[i][j] != INVALID_VALUE) {
                                trace(5, "add_base_value_to_cbias(): pos=%d, bias_cbias=%.3f, base_cbias=%.3f\n",
                                    pos, cbias.cbias[i][pos], temp->cbias[i][j]);
                                cbias.cbias[i][pos] += temp->cbias[i][j];
                                cbias.smode[i][pos] = temp->smode[i][j];
                                cbias.sflag[i][pos] |= 0x01;
                                trace(5, "add_base_value_to_cbias(): pos=%d, bias_cbias=%.3f\n",
                                    pos, cbias.cbias[i][pos]);
                            } else {
                                cbias.smode[i][pos] = temp->smode[i][j];
                                cbias.cbias[i][pos] = INVALID_VALUE;
                                cbias.sflag[i][pos] |= 0x01;
                            }
                        }
                    }
                }
            }
        }
        return &cbias;
    }
    trace(2, "add_base_value_to_cbias(): not found cbias, network=%d, bias_time=%d, base_time=%d\n",
        network, (int)time2gpst(srcbias->time, NULL), (int)time2gpst(basetime, NULL));
    return srcbias;
}

static cssr_bias_bank *add_base_value_to_pbias(cssr_bias_bank *srcbias, int network)
{
    gtime_t basetime = timeadd(srcbias->time, fmod(time2gpst(srcbias->time, NULL), 30.0) * -1.0);
    static cssr_bias_bank pbias;
    cssr_bias_bank *temp;
    int i, j, pos;
    
    if ((temp = get_same_bias_correction(basetime, network)) && (temp->bflag & 0x02)) {
        trace(4, "add_base_value_to_pbias(): bias_tow=%d, network=%d, base_tow=%d, network=%d, flag=0x%02x\n", (int)time2gpst(srcbias->time, NULL),
            srcbias->network, (int)time2gpst(temp->time, NULL), temp->network, temp->bflag);
        memcpy(&pbias, srcbias, sizeof(cssr_bias_bank));
        
        for (i = 0; i < MAXSAT; ++i) {
            if (pbias.prn[i] != 0 && temp->prn[i] != 0) {
                for (j = 0; j < MAXCODE; ++j) {
                    if (temp->smode[i][j] != 0 && (temp->sflag[i][j] & 0x02)) {
                        trace(5, "add_base_value_to_pbias(): j=%2d, prn=%2d, tow=%d, mode=%2d, flag=0x%02x, value=%.3f\n",
                            j, temp->prn[i], (int)time2gpst(temp->time, NULL), temp->smode[i][j],
                            temp->sflag[i][j], temp->pbias[i][j]);
                        if ((pos = get_bias_save_point(&pbias, pbias.prn[i], temp->smode[i][j])) != -1) {
                            if (pbias.pbias[i][pos] != INVALID_VALUE && temp->pbias[i][j] != INVALID_VALUE) {
                                trace(5, "add_base_value_to_pbias(): pos=%d, bias_pbias=%.3f, base_pbias=%.3f\n",
                                    pos, pbias.pbias[i][pos], temp->pbias[i][j]);
                                pbias.pbias[i][pos] += temp->pbias[i][j];
                                pbias.smode[i][pos] = temp->smode[i][j];
                                pbias.sflag[i][pos] |= 0x02;
                                trace(5, "add_base_value_to_pbias(): pos=%d, bias_pbias=%.3f\n",
                                    pos, pbias.pbias[i][pos]);
                            } else {
                                pbias.smode[i][pos] = temp->smode[i][j];
                                pbias.pbias[i][pos] = INVALID_VALUE;
                                pbias.sflag[i][pos] |= 0x02;
                            }
                        }
                    }
                }
            }
        }
        return &pbias;
    }
    trace(2, "add_base_value_to_pbias(): not found pbias, network=%d, bias_time=%d, base_time=%d\n",
        network, (int)time2gpst(srcbias->time, NULL), (int)time2gpst(basetime, NULL));
    return srcbias;
}

static cssr_bias_bank *get_close_network_cbias_correction(gtime_t time, int network, double age)
{
    cssr_bias_bank *cbias = NULL;

    if (network != 0 && !(cbias = get_close_cbias_correction(time, network, age))) {
        if (!(cbias = get_close_cbias_correction(time, 0, age))) {
            return NULL;
        }
    } else if (network == 0) {
        if (!(cbias = get_close_cbias_correction(time, 0, age))) {
            return NULL;
        }
    }
    return cbias;
}

static cssr_bias_bank *get_close_network_pbias_correction(gtime_t time, int network, double age)
{
    cssr_bias_bank *pbias = NULL;
    
    if (network != 0 && !(pbias = get_close_pbias_correction(time, network, age))) {
        if (!(pbias = get_close_pbias_correction(time, 0, age))) {
            return NULL;
        }
    } else if (network == 0) {
        if (!(pbias = get_close_pbias_correction(time, 0, age))) {
            return NULL;
        }
    }
    return pbias;
}

void set_cssr_bank_cbias(gtime_t time, nav_t *nav, int network, int iod)
{
    cssr_bias_bank *bias;
    int i, j, pos;
    
    if ((bias = get_same_bias_correction(time, network)) == NULL) {
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].network = network;
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].bflag = 0x00;
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].time = time;
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].use = TRUE;
        bias = &cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias];
        if (++cssrObject[chidx].NextBias >= BIASBANKNUM) {
            cssrObject[chidx].NextBias = 0;
        }
        memset(bias->smode, 0x00, sizeof(bias->smode));
        memset(bias->sflag, 0x00, sizeof(bias->sflag));
        memset(bias->prn, 0x00, sizeof(bias->prn));
    }

    if (bias->bflag & 0x01) {
        trace(3, "set_cssr_bank_cbias(): clear code bias corrections, network=%d, tow=%d, bflag=0x%02x\n",
            bias->network, (int)time2gpst(bias->time, NULL), bias->bflag);
        for (i = 0; i < MAXSAT; ++i) {
            for (j = 0; j < MAXCODE; ++j) {
                if (bias->sflag[i][j] & 0x01) {
                    bias->sflag[i][j] ^= 0x01;
                    bias->cbias[i][j] = 0.0;
                }
            }
        }
        bias->bflag ^= 0x01;
    }
    
    for (i = 0; i < MAXSAT; ++i) {
        if (nav->ssr[i].update_cb == 1 && nav->ssr[i].nsig > 0) {
            trace(3, "set_cssr_bank_cbias(): tow=%.1f, network=%d, prn=%2d, num=%d\n", time2gpst(time, NULL),
                network, i + 1, nav->ssr[i].nsig);
            bias->prn[i] = i + 1;
            bias->bflag |= 0x01;
            bias->udi[i] = nav->ssr[i].udi[4];
            bias->iod[i] = nav->ssr[i].iod[4];
            
            for (j = 0; j < MAXCODE; ++j) {
                if (nav->ssr[i].smode[j] != 0 && (pos = get_bias_save_point(bias, bias->prn[i], nav->ssr[i].smode[j])) != -1) {
                    bias->cbias[i][pos] = nav->ssr[i].cbias[nav->ssr[i].smode[j]-1];
                    bias->smode[i][pos] = nav->ssr[i].smode[j];
                    bias->sflag[i][pos] |= 0x01;
                    trace(4, "network=%d, pos=%d, mode=%2d, flag=0x%02x, cbias=%.4f\n", network,
                        pos, bias->smode[i][pos], bias->sflag[i][pos], bias->cbias[i][pos]);
                }
            }
        }
    }
    trace(4, "set_cssr_bank_cbias(): next=%d, network=%d, bflag=0x%02x\n", cssrObject[chidx].NextBias, network, bias->bflag);
}

void set_cssr_bank_pbias(gtime_t time, nav_t *nav, int network, int iod)
{
    cssr_bias_bank *bias;
    int i, j, pos;

    if ((bias = get_same_bias_correction(time, network)) == NULL) {
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].network = network;
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].bflag = 0x00;
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].time = time;
        cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias].use = TRUE;
        bias = &cssrObject[chidx].BiasBank[cssrObject[chidx].NextBias];
        if (++cssrObject[chidx].NextBias >= BIASBANKNUM) {
            cssrObject[chidx].NextBias = 0;
        }
        memset(bias->smode, 0x00, sizeof(bias->smode));
        memset(bias->sflag, 0x00, sizeof(bias->sflag));
        memset(bias->prn, 0x00, sizeof(bias->prn));
    }

    if (bias->bflag & 0x02) {
        trace(3, "set_cssr_bank_pbias(): clear phase bias corrections, network=%d, tow=%d, bflag=0x%02x\n",
            bias->network, (int)time2gpst(bias->time, NULL), bias->bflag);
        for (i = 0; i < MAXSAT; ++i) {
            for (j = 0; j < MAXCODE; ++j) {
                if (bias->sflag[i][j] & 0x02) {
                    bias->sflag[i][j] ^= 0x02;
                    bias->pbias[i][j] = 0.0;
                }
            }
        }
        bias->bflag ^= 0x02;
    }
    
    for (i = 0; i < MAXSAT; ++i) {
        if (nav->ssr[i].update_pb == 1 && nav->ssr[i].nsig > 0) {
            trace(3, "set_cssr_bank_pbias(): tow=%.1f, network=%d, prn=%2d, num=%d\n", time2gpst(time, NULL),
                network, i + 1, nav->ssr[i].nsig);
            bias->prn[i] = i + 1;
            bias->bflag |= 0x02;
            bias->udi[i] = nav->ssr[i].udi[5];
            bias->iod[i] = nav->ssr[i].iod[5];
            
            for (j = 0; j < MAXCODE; ++j) {
                if (nav->ssr[i].smode[j] != 0 && (pos = get_bias_save_point(bias, bias->prn[i], nav->ssr[i].smode[j])) != -1) {
                    bias->pbias[i][pos] = nav->ssr[i].pbias[nav->ssr[i].smode[j]-1];
                    bias->smode[i][pos] = nav->ssr[i].smode[j];
                    bias->sflag[i][pos] |= 0x02;
                    trace(4, "network=%d, pos=%d, mode=%2d, flag=0x%02x, pbias=%.4f\n", network,
                        pos, bias->smode[i][pos], bias->sflag[i][pos], bias->pbias[i][pos]);
                }
            }
        }
    }
    trace(4, "set_cssr_bank_pbias(): next=%d, network=%d, bflag=0x%02x\n", cssrObject[chidx].NextBias, network, bias->bflag);
}

static cssr_clock_bank *get_same_clock_correction(gtime_t time, int network)
{
    int i;

    for (i = 0; i < CLKBANKNUM; ++i) {
        if (cssrObject[chidx].ClockBank[i].use == TRUE && cssrObject[chidx].ClockBank[i].network == network &&
            timediff(cssrObject[chidx].ClockBank[i].time, time) == 0.0) {
            return &cssrObject[chidx].ClockBank[i];
        }
    }
    return NULL;
}

static cssr_clock_bank *get_close_clock_correction(gtime_t obstime, gtime_t time, int network, double age)
{
    int search = ((cssrObject[chidx].separation & (1 << (network - 1))) ? network: 0);
    int pos = -1;
    int i;

    for (i = 0; i < CLKBANKNUM; ++i) {
        if (cssrObject[chidx].ClockBank[i].use == TRUE) {
            trace(5, "get_close_clock_correction(): clock=%.1f, network=%d, diff=%.1f\n", time2gpst(cssrObject[chidx].ClockBank[i].time, NULL),
                cssrObject[chidx].ClockBank[i].network, timediff(cssrObject[chidx].ClockBank[i].time, time));
        }
        if (cssrObject[chidx].ClockBank[i].use == TRUE && cssrObject[chidx].ClockBank[i].network == search && timediff(cssrObject[chidx].ClockBank[i].time, time) < age) {
            if (pos != -1 && timediff(cssrObject[chidx].ClockBank[pos].time, cssrObject[chidx].ClockBank[i].time) > 0.0) {
                continue;
            }
            if (timediff(cssrObject[chidx].ClockBank[i].time, obstime) > 0.0) {
                continue;
            }
            pos = i;
        }
    }
    if (pos != -1) {
        trace(4, "get_close_clock_correction(): clock=%.1f, network=%d, diff=%.1f\n", time2gpst(cssrObject[chidx].ClockBank[pos].time, NULL),
            cssrObject[chidx].ClockBank[pos].network, timediff(cssrObject[chidx].ClockBank[pos].time, time));
        return &cssrObject[chidx].ClockBank[pos];
    }
    return NULL;
}

void set_cssr_bank_clock(gtime_t time, nav_t *nav, int network)
{
    cssr_clock_bank *clock;
    int i;

    if ((clock = get_same_clock_correction(time, network)) == NULL) {
        cssrObject[chidx].ClockBank[cssrObject[chidx].NextClock].network = network;
        cssrObject[chidx].ClockBank[cssrObject[chidx].NextClock].time = time;
        cssrObject[chidx].ClockBank[cssrObject[chidx].NextClock].use = TRUE;
        clock = &cssrObject[chidx].ClockBank[cssrObject[chidx].NextClock];
        if (++cssrObject[chidx].NextClock >= CLKBANKNUM) {
            cssrObject[chidx].NextClock = 0;
        }
        memset(clock->prn, 0x00, sizeof(clock->prn));
    }

    for (i = 0; i < MAXSAT; ++i) {
        if (clock->prn[i] == 0) {
            continue;
        }
        clock->prn[i] = 0;
        clock->udi[i] = 0.0;
        clock->iod[i] = 0;
        clock->c0[i] = 0.0;
    }
    
    for (i = 0; i < MAXSAT; ++i) {
        if (nav->ssr[i].update_cc == 1) {
            if (nav->ssr[i].dclk[0] != INVALID_VALUE) {
                trace(3, "set_cssr_bank_clock(): tow=%.1f, network=%d, prn=%2d, c0=%.4f\n", time2gpst(time, NULL),
                    network, i + 1, nav->ssr[i].dclk[0]);
                clock->prn[i] = i + 1;
                clock->udi[i] = nav->ssr[i].udi[1];
                clock->iod[i] = nav->ssr[i].iod[1];
                clock->c0[i] = nav->ssr[i].dclk[0];
            } else {
                trace(3, "set_cssr_bank_clock(): tow=%.1f, network=%d, prn=%2d, c0=#N/A\n", time2gpst(time, NULL),
                    network, i + 1);
                clock->prn[i] = 0;
            }
        }
    }
    trace(4, "set_cssr_bank_clock(): tow=%d, network=%d, next=%d\n",
        (int)time2gpst(time, NULL), network, cssrObject[chidx].NextClock);
}

void set_cssr_latest_trop(gtime_t time, ssrgp_t *ssrg, int network)
{
    int i, j, sat;
    
    for (i = 0; i < ssrg->ngp; ++i) {
        for (j = 0; j < ssrg->nsv[i]; ++j) {
            if (ssrg->stec[i][j] != INVALID_VALUE) {
                sat = ssrg->sat[i][j];
                cssrObject[chidx].LatestTrop.stec0[network-1][i][sat-1] = ssrg->stec0[i][j];
                cssrObject[chidx].LatestTrop.stec[network-1][i][sat-1] = ssrg->stec[i][j];
                cssrObject[chidx].LatestTrop.stectime[network-1][i][sat-1] = time;
                cssrObject[chidx].LatestTrop.prn[network-1][i][sat-1] = sat;
                trace(4, "set_cssr_latest_trop(): network=%d, tow=%.1f, i=%2d, prn=%2d, stec=%.3f, stec0=%.3f\n", network, time2gpst(time, NULL),
                    i, sat, cssrObject[chidx].LatestTrop.stec[network-1][i][sat-1], cssrObject[chidx].LatestTrop.stec0[network-1][i][sat-1]);
            }
        }
        if (ssrg->trop_total[i] != INVALID_VALUE && ssrg->trop_wet[i] != INVALID_VALUE) {
            cssrObject[chidx].LatestTrop.total[network-1][i] = ssrg->trop_total[i];
            cssrObject[chidx].LatestTrop.wet[network-1][i] = ssrg->trop_wet[i];
            cssrObject[chidx].LatestTrop.troptime[network-1][i] = time;
            trace(4, "set_cssr_latest_trop(): network=%d, tow=%.1f, i=%2d, total=%.3f, wet=%.3f\n", network, time2gpst(time, NULL),
                i, cssrObject[chidx].LatestTrop.total[network-1][i], cssrObject[chidx].LatestTrop.wet[network-1][i]);
        }
    }
    cssrObject[chidx].LatestTrop.gridnum[network-1] = ssrg->ngp;
}

int get_cssr_latest_trop(double *total, double *wet, gtime_t time, int network, int index)
{
    if (index < cssrObject[chidx].LatestTrop.gridnum[network-1] && timediff(time, cssrObject[chidx].LatestTrop.troptime[network-1][index]) <= TROPVALIDAGE) {
        trace(4, "get_cssr_latest_trop(): network=%2d, i=%2d, tow=%.1f, diff=%.1f, total=%.3f, wet=%.3f\n", network, index,
            time2gpst(time, NULL), timediff(time, cssrObject[chidx].LatestTrop.troptime[network-1][index]),
            cssrObject[chidx].LatestTrop.total[network-1][index], cssrObject[chidx].LatestTrop.wet[network-1][index]);
        *total = cssrObject[chidx].LatestTrop.total[network-1][index];
        *wet   = cssrObject[chidx].LatestTrop.wet[network-1][index];
        return TRUE;
    }
    trace(4, "get_cssr_latest_trop(): nothing latest trop value, network=%2d, i=%2d, tow=%.1f\n",
        network, index, time2gpst(time, NULL));
    *total = INVALID_VALUE;
    *wet   = INVALID_VALUE;
    return FALSE;
}

int get_cssr_latest_iono(double *iono, gtime_t time, int network, int index, int sat)
{
    if (index < cssrObject[chidx].LatestTrop.gridnum[network-1] && cssrObject[chidx].LatestTrop.prn[network-1][index][sat-1] == sat &&
        timediff(time, cssrObject[chidx].LatestTrop.stectime[network-1][index][sat-1]) <= STECVALIDAGE) {
        if (cssrObject[chidx].LatestTrop.stec[network-1][index][sat-1] == INVALID_VALUE) {
            *iono = 40.3E16 / (FREQ1 * FREQ2) * cssrObject[chidx].LatestTrop.stec0[network-1][index][sat-1];
        } else {
            *iono = 40.3E16 / (FREQ1 * FREQ2) * cssrObject[chidx].LatestTrop.stec[network-1][index][sat-1];
        }
        trace(4, "get_cssr_latest_iono(): network=%2d, i=%2d, sat=%2d, tow=%.1f, diff=%.1f, iono=%.3f\n", network, index,
            sat, time2gpst(time, NULL), timediff(time, cssrObject[chidx].LatestTrop.stectime[network-1][index][sat-1]), *iono);
        return TRUE;
    }
    *iono = INVALID_VALUE;
    return FALSE;
}

static cssr_trop_bank *get_same_trop_correction(gtime_t time)
{
    int i;

    for (i = 0; i < TROPBANKNUM; ++i) {
        if (cssrObject[chidx].TropBank[i].use == TRUE && timediff(cssrObject[chidx].TropBank[i].time, time) == 0.0) {
            return &cssrObject[chidx].TropBank[i];
        }
    }
    return NULL;
}

static cssr_trop_bank *get_close_trop_correction(gtime_t obstime, gtime_t time, int network, double age, int tropless)
{
    int pos = -1;
    int i;

    trace(4, "get_close_trop_correction(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);

    if (network >= 1 && network < CSSR_MAX_NETWORK) {
        for (i = 0; i < TROPBANKNUM; ++i) {
            if (cssrObject[chidx].TropBank[i].use == TRUE) {
                trace(5, "get_close_trop_correction(): trop=%.1f, diff=%.1f, ngp=%d\n", time2gpst(cssrObject[chidx].TropBank[i].time, NULL),
                    timediff(cssrObject[chidx].TropBank[i].time, time), cssrObject[chidx].TropBank[i].gridnum[network-1]);
            }
            if (cssrObject[chidx].TropBank[i].use == TRUE && fabs(timediff(time, cssrObject[chidx].TropBank[i].time)) <= age) {
                if (pos != -1 && timediff(cssrObject[chidx].TropBank[pos].time, cssrObject[chidx].TropBank[i].time) > 0.0) {
                    continue;
                }
                if (tropless && timediff(obstime, cssrObject[chidx].TropBank[i].time) < 0.0) {
                    continue;
                }
                if (cssrObject[chidx].TropBank[i].gridnum[network-1] > 0) {
                    pos = i;
                }
                trace((network <= 12 ? 3: 5), "get_close_trop_correction(): network=%2d, obs=%.1f, trop=%d, grids=%2d\n",
                    network, time2gpst(obstime, NULL), (int)time2gpst(cssrObject[chidx].TropBank[i].time, NULL),
                    cssrObject[chidx].TropBank[i].gridnum[network-1]);
            }
        }
        if (pos != -1) {
            trace(4, "get_close_trop_correction(): trop=%.1f, diff=%.1f, ngp=%d\n", time2gpst(cssrObject[chidx].TropBank[pos].time, NULL),
                timediff(cssrObject[chidx].TropBank[pos].time, time), cssrObject[chidx].TropBank[pos].gridnum[network-1]);
            return &cssrObject[chidx].TropBank[pos];
        }
    }
    return NULL;
}

void set_cssr_bank_trop(gtime_t time, ssrgp_t *ssrg, int network)
{
    cssr_trop_bank *trop;
    int i, j;

    if ((trop = get_same_trop_correction(time)) == NULL) {
        cssrObject[chidx].TropBank[cssrObject[chidx].NextTrop].time = time;
        cssrObject[chidx].TropBank[cssrObject[chidx].NextTrop].use = TRUE;
        trop = &cssrObject[chidx].TropBank[cssrObject[chidx].NextTrop];
        if (++cssrObject[chidx].NextTrop >= TROPBANKNUM) {
            cssrObject[chidx].NextTrop = 0;
        }
        memset(trop->gridnum, 0x00, sizeof(trop->gridnum));
        memset(trop->satnum, 0x00, sizeof(trop->satnum));
    }

    trop->gridnum[network-1] = ssrg->ngp;
    trace(4, "set_cssr_bank_trop(): network=%d, tow=%.1f, gridnum=%d\n", network,
        time2gpst(trop->time, NULL), trop->gridnum[network-1]);
    for (i = 0; i < ssrg->ngp; ++i) {
        trop->gridpos[network-1][i][0] = ssrg->gp[i].pos[0];
        trop->gridpos[network-1][i][1] = ssrg->gp[i].pos[1];
        trop->gridpos[network-1][i][2] = ssrg->gp[i].pos[2];
        trop->total[network-1][i] = ssrg->trop_total[i];
        trop->wet[network-1][i] = ssrg->trop_wet[i];
        trop->satnum[network-1][i] = ssrg->nsv[i];
        if (trop->total[network-1][i] == INVALID_VALUE || trop->wet[network-1][i] == INVALID_VALUE) {
            get_cssr_latest_trop(&trop->total[network-1][i], &trop->wet[network-1][i], time, network, i);
        }
        trace(4, "set_cssr_bank_trop(): i=%d, lat=%.4f, lon=%.4f, alt=%.1f, total=%.3f, wet=%.3f, satnum=%d\n", i, trop->gridpos[network-1][i][0]*R2D, trop->gridpos[network-1][i][1]*R2D,
                trop->gridpos[network-1][i][2], trop->total[network-1][i], trop->wet[network-1][i], trop->satnum[network-1][i]);
        for (j = 0; j < trop->satnum[network-1][i]; ++j) {
            if (ssrg->stec[i][j] == INVALID_VALUE) {
                get_cssr_latest_iono(&trop->iono[network-1][i][j], time, network, i, ssrg->sat[i][j]);
            } else {
                trop->iono[network-1][i][j] = 40.3E16 / (FREQ1 * FREQ2) * ssrg->stec[i][j];
            }
            trop->prn[network-1][i][j] = ssrg->sat[i][j];
            trace(4, "set_cssr_bank_trop(): prn=%d, iono=%.3f\n", trop->prn[network-1][i][j], trop->iono[network-1][i][j]);
        }
    }
    trace(4, "set_cssr_bank_trop(): next=%d\n", cssrObject[chidx].NextTrop);
}

void set_grid_data(double *pos, int index, int chindex)
{
    CurrentCSSR[chindex].stec[index].pos[0] = pos[0] * R2D;
    CurrentCSSR[chindex].stec[index].pos[1] = pos[1] * R2D;
    CurrentCSSR[chindex].stec[index].pos[2] = pos[2];
    CurrentCSSR[chindex].stec[index].n = 0;

    CurrentCSSR[chindex].zwd[index].pos[0] = (float)(pos[0] * R2D);
    CurrentCSSR[chindex].zwd[index].pos[1] = (float)(pos[1] * R2D);
    CurrentCSSR[chindex].zwd[index].pos[2] = (float)pos[2];
    CurrentCSSR[chindex].zwd[index].n = 0;
}

void init_grid_index(cssr_bank *bank)
{
    int i;

    for (i = 0; i < RTCM_SSR_MAX_GP; ++i) {
        bank->stec[i].data = &bank->stecdata[i][0];
        bank->stec[i].nmax = MAXSAT;
        bank->zwd[i].data = &bank->zwddata[i];
        bank->zwd[i].nmax = 1;
    }
}

extern void init_fastfix_flag(void)
{
    int i;
    
    for (i = 0; i < CSSR_MAX_NETWORK - 1; ++i) {
        cssrObject[chidx].fastfix[i+1] = TRUE;
    }
}

extern void check_cssr_grid_status(gtime_t time)
{
    int i, j, k, valid, network, preliminary;
    double tow = time2gpst(time, NULL);
    cssr_orbit_bank *orbit = NULL;
    cssr_trop_bank *trop = NULL;
    char dbgstr[144];
    
    for (i = 0; i < CSSR_MAX_NETWORK - 1; ++i) {
        for (j = 0, preliminary = i + 1; j < 2; ++j) {
            if ((orbit = get_close_orbit_correction(time, (j == 0 ? preliminary: 0), 180.0)) != NULL) {
                network = (orbit->network ? orbit->network: preliminary);
                break;
            }
        }
        if (orbit && !(trop = get_close_trop_correction(time, orbit->time, network, 30.0, TRUE)) &&
            cssrObject[chidx].fastfix[network] == TRUE) {
            trop = get_close_trop_correction(time, orbit->time, network, 30.0, FALSE);
        }
        
        if (!orbit || !trop) {
            sprintf(dbgstr, "check_cssr_grid_status(): correction is not found, network=%02d",
                preliminary);
            sprintf(&dbgstr[strlen(dbgstr)], ", obs=%.1f", time2gpst(time, NULL));
            if (orbit) {
                sprintf(&dbgstr[strlen(dbgstr)], ", orbit=%d", (int)time2gpst(orbit->time, NULL));
            } else {
                strcat(dbgstr, ", orbit=NULL");
            }
            if (trop) {
                sprintf(&dbgstr[strlen(dbgstr)], ", trop=%d\n", (int)time2gpst(trop->time, NULL));
            } else {
                strcat(dbgstr, ", trop=NULL\n");
            }
            trace((preliminary <= 12 ? 2: 5), dbgstr);
            for (j = 0; j < RTCM_SSR_MAX_GP; ++j) {
                grid_stat[chidx][preliminary][j] = FALSE;
            }
            continue;
        }
        if (cssrObject[chidx].fastfix[network] == TRUE && timediff(time, trop->time) >= 0.0) {
            cssrObject[chidx].fastfix[network] = FALSE;
        }
        
        trace(4, "check_cssr_grid_status(): network=%2d, fastfix=%d, obs=%.1f, trop=%d\n",
            network, cssrObject[chidx].fastfix[network], tow, (int)time2gpst(trop->time, NULL));
        for (j = 0; j < trop->gridnum[i]; ++j) {
            if (trop->total[i][j] == INVALID_VALUE || trop->wet[i][j] == INVALID_VALUE) {
                trace(3, "%2d-%2d, total=%.4f, wet=%.4f\n", i, j, trop->total[i][j],
                    trop->wet[i][j]);
                grid_stat[chidx][network][j] = FALSE;
                continue;
            }
            for (k = 0, valid = trop->satnum[i][j]; k < trop->satnum[i][j]; ++k) {
                if (trop->iono[i][j][k] == INVALID_VALUE) {
                    valid -= 1;
                    continue;
                }
            }
            if (valid < 8) {
                trace(3, "%2d-%2d, valid=%2d\n", i, j, valid);
                grid_stat[chidx][network][j] = FALSE;
                continue;
            }
            grid_stat[chidx][network][j] = TRUE;
        }
    }
}

static int sub_get_close_cssr(gtime_t time, int network, cssr_orbit_bank **orbit, cssr_clock_bank **clock, cssr_bias_bank **cbias, cssr_bias_bank **pbias, cssr_trop_bank **trop, int *flag)
{
    double tow;
    
    trace(4, "sub_get_close_cssr(): tow=%.1f, network=%d\n", (tow = time2gpst(time, NULL)), network);
    
    if (!(*orbit = get_close_orbit_correction(time, network, 180.0)) || !(*pbias = get_close_network_pbias_correction((*orbit)->time, network, 0.0))) {
        trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, orbit or pbias correction is not found.\n", network, tow);
        return FALSE;
    }
    if (!(*trop = get_close_trop_correction(time, (*orbit)->time, network, 30.0, (cssrObject[chidx].fastfix[network] ? FALSE: TRUE)))) {
        trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, orbit=%d, trop correction is not found.\n",
            network, tow, (int)time2gpst((*orbit)->time, NULL));
        return FALSE;
    }
    
    if (timediff((*trop)->time, (*pbias)->time) < 0.0 || timediff((*trop)->time, (*pbias)->time) >= 30.0) {
        if (!(*pbias = get_close_network_pbias_correction((*trop)->time, network, 0.0))) {
            trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, trop=%d, pbias correction is not found.\n",
                network, tow, (int)time2gpst((*trop)->time, NULL));
            return FALSE;
        }
    }
    if (!(*cbias = get_close_network_cbias_correction((*pbias)->time, network, 0.0))) {
        trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, pbias=%d, cbias correction is not found.\n",
            network, tow, (int)time2gpst((*pbias)->time, NULL));
        return FALSE;
    }
    
    if (!(*clock = get_close_clock_correction(time, (*orbit)->time, network, 30.0))) {
        trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, clock correction is not found.\n",
            network, tow);
        return FALSE;
    }
    
    if (timediff((*pbias)->time, (*cbias)->time) < -30.0 || timediff((*pbias)->time, (*cbias)->time) >  0.0) {
        trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, cbias=%d, pbias=%d, time is not much.\n", network,
            tow, (int)time2gpst((*cbias)->time, NULL), (int)time2gpst((*pbias)->time, NULL));
        return FALSE;
    }
    if (timediff((*trop)->time, (*orbit)->time) <  -30.0 || timediff((*trop)->time, (*orbit)->time) >= 30.0) {
        trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, orbit=%d, trop=%d, time is not much.\n", network,
            tow, (int)time2gpst((*orbit)->time, NULL), (int)time2gpst((*trop)->time, NULL));
        return FALSE;
    }
    if (timediff((*trop)->time, (*pbias)->time) <    0.0 || timediff((*trop)->time, (*pbias)->time) >= 30.0) {
        trace(2, "sub_get_close_cssr(): network=%d, obs=%.1f, pbias=%d, trop=%d, time is not much.\n", network,
            tow, (int)time2gpst((*pbias)->time, NULL), (int)time2gpst((*trop)->time, NULL));
        return FALSE;
    }
    
    *cbias = add_base_value_to_cbias(*cbias, ((*cbias)->network == 0 ? network: 0));
    *pbias = add_base_value_to_pbias(*pbias, ((*pbias)->network == 0 ? network: 0));
    *flag = (timediff((*cbias)->time, (*pbias)->time) == 0.0 ? TRUE: FALSE);
    trace(5, "sub_get_close_cssr(): obs=%.1f, flag=%d\n", tow, *flag);
    return TRUE;
}

extern int get_close_cssr(gtime_t time, int network, int l6mrg)
{
    static cssr_bank saveCSSR[L6_CH_NUM];
    cssr_orbit_bank *orbit;
    cssr_clock_bank *clock;
    cssr_bias_bank *cbias;
    cssr_bias_bank *pbias;
    cssr_trop_bank *trop;
    int i, j, l, ch;
    int sis, ret[SSR_CH_NUM]={0}, r=0;

    for (ch=0;ch<(l6mrg?SSR_CH_NUM:1);ch++) {
        set_cssr_ch_idx(ch);
        if (!cssrObject[ch].use) {
            break;
        }
        if ((ret[ch]=sub_get_close_cssr(time, network, &orbit, &clock, &cbias, &pbias, &trop, &sis)) == FALSE) {
            if (ch==(l6mrg?SSR_CH_NUM:1)-1) {
                for (l=0;l<SSR_CH_NUM;l++) { r+=ret[l]; }
                if (r==0) { return FALSE; } else { return TRUE; }
            } else {
                continue;
            }
        }

        trace(3, "update cssr: facility=%d, network=%d, obs=%.1f, orbit=%.1f, clock=%.1f, cbias=%.1f, pbias=%.1f, trop=%.1f, gridnum=%d\n", cssrObject[ch].Facility + 1,
            network, time2gpst(time, NULL), time2gpst(orbit->time, NULL), time2gpst(clock->time, NULL), time2gpst(cbias->time, NULL),
            time2gpst(pbias->time, NULL), time2gpst(trop->time, NULL), trop->gridnum[network - 1]);
        memcpy(&saveCSSR[ch], &CurrentCSSR[ch], sizeof(CurrentCSSR[ch]));
        memset(&CurrentCSSR[ch], 0x00, sizeof(CurrentCSSR[ch]));
        init_grid_index(&CurrentCSSR[ch]);
        init_grid_index(&saveCSSR[ch]);

        for (i = 0; i < MAXSAT; ++i) {
            if (cbias->prn[i] != 0) {
                for (j = 0; j < MAXCODE; ++j) {
                    if (cbias->smode[i][j] != 0 && (cbias->sflag[i][j] & 0x01) == 0x01) {
                        trace(4, "sub_get_close_cssr(): cbias, network=%d, j=%d, prn=%d, mode=%d, cbias=%.6f\n", network, j,
                            cbias->prn[i], cbias->smode[i][j], cbias->cbias[i][j]);
                        CurrentCSSR[ch].cbias[i][j] = cbias->cbias[i][j];
                        CurrentCSSR[ch].smode[i][j] = cbias->smode[i][j];
                    }
                }
            }
            if (pbias->prn[i] != 0) {
                for (j = 0; j < MAXCODE; ++j) {
                    if (pbias->smode[i][j] != 0 && (pbias->sflag[i][j] & 0x02) == 0x02) {
                        trace(4, "sub_get_close_cssr(): pbias, network=%d, j=%d, prn=%d, mode=%d, pbias=%.6f\n", network, j,
                            pbias->prn[i], pbias->smode[i][j], pbias->pbias[i][j]);
                        CurrentCSSR[ch].pbias[i][j] = pbias->pbias[i][j];
                        CurrentCSSR[ch].smode[i][j] = pbias->smode[i][j];
                    }
                }
            }
            CurrentCSSR[ch].time[i][4] = cbias->time;
            CurrentCSSR[ch].prn[i][4] = cbias->prn[i];
            CurrentCSSR[ch].udi[i][4] = cbias->udi[i];
            CurrentCSSR[ch].iod[i][4] = cbias->iod[i];
            CurrentCSSR[ch].time[i][5] = pbias->time;
            CurrentCSSR[ch].prn[i][5] = pbias->prn[i];
            CurrentCSSR[ch].udi[i][5] = pbias->udi[i];
            CurrentCSSR[ch].iod[i][5] = pbias->iod[i];
            CurrentCSSR[ch].flag[i] |= 0x04;
        }

        for (i = 0; i < MAXSAT; ++i) {
            if (orbit->prn[i] != 0) {
                trace(4, "sub_get_close_cssr(): orbit, prn=%d, iode=%d, deph0=%f, deph1=%f, deph2=%f\n", orbit->prn[i],
                    orbit->iode[i], orbit->deph0[i], orbit->deph1[i], orbit->deph2[i]);
                CurrentCSSR[ch].time[i][0] = orbit->time;
                CurrentCSSR[ch].prn[i][0] = orbit->prn[i];
                CurrentCSSR[ch].udi[i][0] = orbit->udi[i];
                CurrentCSSR[ch].iod[i][0] = orbit->iod[i];
                CurrentCSSR[ch].iode[i] = orbit->iode[i];
                CurrentCSSR[ch].deph0[i] = orbit->deph0[i];
                CurrentCSSR[ch].deph1[i] = orbit->deph1[i];
                CurrentCSSR[ch].deph2[i] = orbit->deph2[i];
                CurrentCSSR[ch].flag[i] |= 0x01;
            }
        }

        for (i = 0, CurrentCSSR[ch].gridnum = trop->gridnum[network - 1]; i < CurrentCSSR[ch].gridnum; ++i) {
            for (j = 0, set_grid_data(&trop->gridpos[network - 1][i][0], i, ch); j < trop->satnum[network - 1][i]; ++j) {
                if (trop->iono[network - 1][i][j] != INVALID_VALUE) {
                    add_data_stec(&CurrentCSSR[ch].stec[i], trop->time, trop->prn[network - 1][i][j],
                        0, trop->iono[network - 1][i][j], 0.0, 0.0, 0.0);
                }
            }
            add_data_trop(&CurrentCSSR[ch].zwd[i], trop->time, trop->wet[network - 1][i],
                trop->total[network - 1][i], 9999, 0, 1);
        }

        for (i = 0; i < MAXSAT; ++i) {
            if (clock->prn[i] != 0) {
                trace(4, "sub_get_close_cssr(): clock, prn=%d, c0=%f\n", clock->prn[i], clock->c0[i]);
                CurrentCSSR[ch].time[i][1] = clock->time;
                CurrentCSSR[ch].prn[i][1] = clock->prn[i];
                CurrentCSSR[ch].udi[i][1] = clock->udi[i];
                CurrentCSSR[ch].iod[i][1] = clock->iod[i];
                CurrentCSSR[ch].c0[i] = clock->c0[i];
                CurrentCSSR[ch].flag[i] |= 0x02;
            }
        }

        CurrentCSSR[ch].orbit_time = orbit->time;
        CurrentCSSR[ch].clock_time = clock->time;
        CurrentCSSR[ch].bias_time = cbias->time;
        CurrentCSSR[ch].trop_time = trop->time;
        CurrentCSSR[ch].facility = cssrObject[ch].Facility;
        CurrentCSSR[ch].update_time = time;
        CurrentCSSR[ch].network = network;
        CurrentCSSR[ch].use = TRUE;
    }
    return TRUE;
}

extern void update_global_cssr(ssr_t *ssr, int sat, int ch)
{
    int i;

    if (CurrentCSSR[ch].prn[sat-1][4] != 0) {
            for (i = ssr->nsig = 0; i < MAXCODE; ++i) {
            ssr->discontinuity[i] = 0;
            ssr->pbias[i] = 0.0;
            ssr->cbias[i] = 0.0;
            ssr->smode[i] = 0;
        }
        for (i = 0; i < MAXCODE; ++i) {
            if (CurrentCSSR[ch].smode[sat-1][i] != 0) {
                trace(3, "update_global_cssr(): bias correction, network=%d, tow=%.1f, sat=%d, mode=%d, cbias=%.6f, pbias=%.6f\n", CurrentCSSR[ch].network,
                        time2gpst(CurrentCSSR[ch].time[sat-1][4], NULL), sat, CurrentCSSR[ch].smode[sat-1][i],
                        CurrentCSSR[ch].cbias[sat-1][i], CurrentCSSR[ch].pbias[sat-1][i]);
                ssr->cbias[CurrentCSSR[ch].smode[sat-1][i]-1] = CurrentCSSR[ch].cbias[sat-1][i];
                ssr->pbias[CurrentCSSR[ch].smode[sat-1][i]-1] = CurrentCSSR[ch].pbias[sat-1][i];
                ssr->smode[i] = CurrentCSSR[ch].smode[sat-1][i];
                ++ssr->nsig;
            }
        }
        ssr->udi[4] = CurrentCSSR[ch].udi[sat-1][4];
        ssr->udi[5] = CurrentCSSR[ch].udi[sat-1][5];
        ssr->iod[4] = CurrentCSSR[ch].iod[sat-1][4];
        ssr->iod[5] = CurrentCSSR[ch].iod[sat-1][5];
        ssr->t0[4] = CurrentCSSR[ch].time[sat-1][4];
        ssr->t0[5] = CurrentCSSR[ch].time[sat-1][5];
    } else {
        memset(&ssr->t0[4], 0x00, sizeof(ssr->t0[4]));
        memset(&ssr->t0[5], 0x00, sizeof(ssr->t0[5]));
    }
    
    if (CurrentCSSR[ch].prn[sat-1][0] != 0) {
        trace(3, "update_global_cssr(): orbit correction, network=%d, tow=%.1f, sat=%d, iode=%d, deph0=%.6f\n", CurrentCSSR[ch].network,
                time2gpst(CurrentCSSR[ch].time[sat-1][0], NULL), sat, CurrentCSSR[ch].iode[sat-1], CurrentCSSR[ch].deph0[sat-1]);
        ssr->t0[0] = CurrentCSSR[ch].time[sat-1][0];
        ssr->udi[0] = CurrentCSSR[ch].udi[sat-1][0];
        ssr->iod[0] = CurrentCSSR[ch].iod[sat-1][0];
        ssr->iode = CurrentCSSR[ch].iode[sat-1];
        ssr->deph[0] = CurrentCSSR[ch].deph0[sat-1];
        ssr->deph[1] = CurrentCSSR[ch].deph1[sat-1];
        ssr->deph[2] = CurrentCSSR[ch].deph2[sat-1];
        ssr->ddeph[0] = 0.0;
        ssr->ddeph[1] = 0.0;
        ssr->ddeph[2] = 0.0;
    } else {
        memset(&ssr->t0[0], 0x00, sizeof(ssr->t0[0]));
    }
    
    if (CurrentCSSR[ch].prn[sat-1][1] != 0) {
        trace(3, "update_global_cssr(): clock correction, network=%d, tow=%.1f, sat=%d, dclk=%.6f\n", CurrentCSSR[ch].network,
                time2gpst(CurrentCSSR[ch].time[sat-1][1], NULL), sat, CurrentCSSR[ch].c0[sat-1]);
        ssr->t0[1] = CurrentCSSR[ch].time[sat-1][1];
        ssr->udi[1] = CurrentCSSR[ch].udi[sat-1][1];
        ssr->iod[1] = CurrentCSSR[ch].iod[sat-1][1];
        ssr->dclk[0] = CurrentCSSR[ch].c0[sat-1];
        ssr->dclk[1] = ssr->dclk[2] = 0.0;
        ssr->update = 1;
    } else {
        memset(&ssr->t0[1], 0x00, sizeof(ssr->t0[1]));
    }
}

extern void update_local_cssr(nav_t *nav, int l6mrg)
{
    int ch;
    for (ch=0;ch<(l6mrg?SSR_CH_NUM:1);ch++) {
        nav->stec_ch[ch] = &CurrentCSSR[ch].stec[0];
        nav->zwd_ch[ch]  = &CurrentCSSR[ch].zwd[0];
    }        
}

extern void check_cssr_facility(nav_t *nav, int network, int l6mrg)
{
    static int savefacility[L6_CH_NUM] = {-1, -1};
    static int savenetwork[L6_CH_NUM] = {-1, -1};
    static int savedelivery[L6_CH_NUM] = {-1, -1};
    int ch;

    for (ch=0;ch<(l6mrg?SSR_CH_NUM:1);ch++) {
        if (savefacility[ch] != CurrentCSSR[ch].facility) {
            if (savefacility[ch] != -1) {
                trace(1, "CSSR facility changed[ch:%d] %d(%d) ---> %d(%d)\n", ch, savefacility[ch] + 1, savedelivery[chidx], CurrentCSSR[ch].facility + 1, l6delivery[chidx]);
                nav->filreset = TRUE;
            } else {
                trace(1, "CSSR facility changed[ch:%d]        ---> %d(%d)\n", ch, CurrentCSSR[ch].facility + 1, l6delivery[chidx]);
            }
        }
        if (savenetwork[ch] != -1 && savenetwork[ch] != network) {
            trace(1, "CSSR network changed[ch:%d], %d ---> %d\n", ch, savenetwork[ch], network);
            if (nav->filreset != TRUE) {
                nav->ionoreset = TRUE;
            }
        }
        savedelivery[chidx] = l6delivery[chidx];
        savefacility[ch] = CurrentCSSR[ch].facility;
        savenetwork[ch] = network;
    }
}

extern int get_current_cssr_facility(int idx)
{
    return(CurrentCSSR[idx].facility + 1);
}

gtime_t get_backup_cssr_time(int idx)
{
    gtime_t time = BackupCSSR[idx].clock_time;

    if (timediff(BackupCSSR[idx].orbit_time, time) > 0.0) {
        time = BackupCSSR[idx].orbit_time;
    }
    if (timediff(BackupCSSR[idx].bias_time, time) > 0.0) {
        time = BackupCSSR[idx].bias_time;
    }
    return time;
}

extern int is_valid_cssr_backup(gtime_t obstime, int l6mrg){
    int ch;

    for (ch=0;ch<(l6mrg?SSR_CH_NUM:1);ch++) {
        if (timediff(obstime, get_backup_cssr_time(ch)) <= 180.0) {
            return TRUE;
        }    
    }
    return FALSE;
}

extern void backup_current_cssr(grid_t *grid, int l6mrg)
{
    int ch;

    for (ch=0;ch<(l6mrg?SSR_CH_NUM:1);ch++) {
        memcpy(&BackupCSSR[ch], &CurrentCSSR[ch], sizeof(cssr_bank));
        memcpy(&BackupGrid[ch], grid, sizeof(grid_t));
        init_grid_index(&BackupCSSR[ch]);
    }
}

extern void restore_current_cssr(gtime_t time, grid_t *grid, int l6mrg)
{
    cssr_clock_bank *clock;
    int i,ch;
    
    for (ch=0;ch<(l6mrg?SSR_CH_NUM:1);ch++) {
        set_cssr_ch_idx(ch);
        if (!cssrObject[ch].use) {
            break;
        }
        if (timediff(time, get_backup_cssr_time(ch)) > 180.0) { continue; }

        memcpy(&CurrentCSSR[ch], &BackupCSSR[ch], sizeof(cssr_bank));
        memcpy(grid, &BackupGrid[ch], sizeof(grid_t));
        init_grid_index(&CurrentCSSR[ch]);

        if ((clock = get_close_clock_correction(time, CurrentCSSR[ch].orbit_time, grid->network, 30.0))) {
            if (timediff(CurrentCSSR[ch].clock_time, clock->time) != 0.0) {
                trace(2, "restore_current_cssr(): network=%2d, diff=%.1f, orbit=%d, clock=%d\n",
                    grid->network, timediff(clock->time, CurrentCSSR[ch].orbit_time),
                    (int)time2gpst(CurrentCSSR[ch].orbit_time, NULL),
                    (int)time2gpst(clock->time, NULL));
                CurrentCSSR[ch].clock_time = clock->time;

                for (i = 0; i < MAXSAT; ++i) {
                    if (clock->prn[i] != 0) {
                        trace(4, "restore_current_cssr(): clock, prn=%2d, c0=%.3f\n",
                            clock->prn[i], clock->c0[i]);
                        CurrentCSSR[ch].time[i][1] = clock->time;
                        CurrentCSSR[ch].prn[i][1] = clock->prn[i];
                        CurrentCSSR[ch].udi[i][1] = clock->udi[i];
                        CurrentCSSR[ch].iod[i][1] = clock->iod[i];
                        CurrentCSSR[ch].c0[i] = clock->c0[i];
                        CurrentCSSR[ch].flag[i] |= 0x02;
                    }
                    else {
                        CurrentCSSR[ch].flag[i] ^= 0x02;
                        CurrentCSSR[ch].prn[i][1] = 0;
                    }
                }
            }
        }
    }
}

extern void clear_current_cssr(void)
{
    int i;
    for (i = 0; i < L6_CH_NUM; i++) {
        memset(&BackupCSSR[i], 0x00, sizeof(BackupCSSR[i]));
        init_grid_index(&BackupCSSR[i]);
        memset(&CurrentCSSR[i], 0x00, sizeof(cssr_bank));
        init_grid_index(&CurrentCSSR[i]);
    }
}

static void output_cssr_head(rtcm_t *rtcm, cssr_t *cssr, int sync, int tow, int udi, int iod, int ngnss,  FILE *fp)
{
    if (fp == NULL) return;

    fprintf(fp, "%d, %d, %d, %d, %d, %d", rtcm->ctype, rtcm->subtype, tow, udi, sync, iod);
    if (rtcm->subtype == CSSR_TYPE_MASK) {
        fprintf(fp, ", %d", ngnss);
    }
}

/* decode cssr header ---------------------------------------------------------*/
static int decode_cssr_head(rtcm_t *rtcm, cssr_t *cssr, int *sync, int *tow,
                           int *iod, int *iod_sv, double *udint, int *ngnss, int i0, int header, FILE *fp)
{
    int i=i0,udi;

    if (rtcm->subtype ==  CSSR_TYPE_MASK) {
        *tow = getbitu(rtcm->buff,i,20); i+=20; /* gps epoch time */
    } else {
        *tow = rtcm->tow0 + getbitu(rtcm->buff,i,12); i+=12; /* gps epoch time (hourly) */
    }

    trace(4,"decode_cssr_head: subtype=%d epoch=%4d\n", rtcm->subtype, *tow);

    udi = getbitu(rtcm->buff,i, 4); i+= 4; /* update interval */
    *sync = getbitu(rtcm->buff,i, 1); i+= 1; /* multiple message indicator */
    *udint = ssrudint[udi];
    *iod = getbitu(rtcm->buff,i, 4); i+= 4; /* iod ssr */

    if (rtcm->subtype == CSSR_TYPE_MASK) {
        cssr->iod = *iod;
        *ngnss = getbitu(rtcm->buff,i, 4); i+= 4; /* number of gnss */
    }
    output_cssr_head(rtcm, cssr, *sync, *tow, udi, *iod, *ngnss, fp);

    return i;
}

static void output_cssr_mask(rtcm_t *rtcm, cssr_t *cssr, int ngnss, int *id_, FILE *fp) {
    int nsat_g=0, ncell=0, prn, sat[CSSR_MAX_SV];
    int j, k, id;

    if (fp == NULL) return;

    trace(4, "output_cssr_mask: ngnss=%d\n", ngnss);
    for (k=0;k<ngnss;k++) {
        if (k != 0) {
            fprintf(fp, ",,,,,,, ");
        } else {
            fprintf(fp, ", ");
        }
        id = id_[k];    
        fprintf(fp, "%d, ", id);
        fprintf(fp, "0x%02x",  (uint32_t)(cssr->svmask[id]>>32));
        fprintf(fp, "%08lx, ", (uint64_t)cssr->svmask[id] & 0x00000000ffffffff);
        fprintf(fp, "0x%04x, ", cssr->sigmask[id]);
        fprintf(fp, "%d", cssr->cmi[id]);
        nsat_g = svmask2nsatlist(cssr->svmask[id], id, sat);

        if (cssr->cmi[id]) { /* cell-mask is included */
            for (j = 0; j < nsat_g; j++) {
                satsys(sat[j], &prn);
                if (j != 0) {
                    fprintf(fp, ",,,,,,,,,,, %d, 0x%04x\n", prn, cssr->cellmask[ncell]);
                } else {
                    fprintf(fp, ", %d, 0x%04x\n", prn, cssr->cellmask[ncell]);
                }
                ncell++;
            }
        } else {
            fprintf(fp, "\n");
        }
    }
    if (ngnss <= 0) {
        fprintf(fp, "\n");
    }
}

/* decode mask message */
static int decode_cssr_mask(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, l, sync, tow, ngnss, iod, nsat_g=0, id, nsig, ncell=0, prn, sat[CSSR_MAX_SV];
    int id_[CSSR_MAX_GNSS];
    double udint;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    rtcm->tow0 = floor(tow/3600.0)*3600.0;
    check_week_ref(rtcm, tow, ref_mask);
    rtcm->time = gpst2time(rtcm->week_ref[ref_mask], tow);
    for (j = 0; j < CSSR_MAX_GNSS; j++) {
        cssr->cmi[j] = 0;
        cssr->svmask[j] = 0;
        cssr->sigmask[j] = 0;
    }
    for (j = 0; j < CSSR_MAX_SV; j++) {
        cssr->cellmask[j] = 0;
    }

    trace(2,"decode_cssr_mask: facility=%d tow=%d iod=%d\n", l6facility[chidx]+1, tow, cssr->iod);

    for (k=0;k<ngnss;k++) {
        id_[k] = id = getbitu(rtcm->buff,i, 4); i+= 4; /* gnss id */
        cssr->svmask[id] = (uint64_t)getbitu(rtcm->buff,i, 8)<<32; i+= 8; /* sv mask */
        cssr->svmask[id] |= getbitu(rtcm->buff,i, 32); i+= 32; /* sv mask */
        cssr->sigmask[id] = getbitu(rtcm->buff,i, 16); i+= 16; /* signal mask */
        cssr->cmi[id] = getbitu(rtcm->buff,i, 1); i++; /* cell mask availability */

        nsig = sigmask2nsig(cssr->sigmask[id]);
        nsat_g = svmask2nsatlist(cssr->svmask[id], id, sat);

        if (cssr->cmi[id]) { /* cell-mask is included */
            for (j = 0; j < nsat_g; j++) {
                cssr->cellmask[ncell] = getbitu(rtcm->buff, i, nsig); i += nsig;
                satsys(sat[j], &prn);
                ncell++;
            }
        } else {
            for (j = 0; j < nsat_g; j++) {
                for (l = 0; l < nsig; l++) {
                    cssr->cellmask[ncell] |= ((uint16_t)1<<(nsig-1-l));
                }
                ++ncell;
            }
        }
    }
    output_cssr_mask(rtcm, cssr, ngnss, id_, fp);
    cssr->l6delivery = l6delivery[chidx];
    cssr->l6facility = l6facility[chidx];
    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is enough to decode the mask message */
static int check_bit_width_mask(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int k,ngnss=0,cmi=0,nsig,nsat;
    uint16_t sigmask;
    uint64_t svmask;

    if (i0+49>rtcm->havebit) return FALSE;

    ngnss = getbitu(rtcm->buff, i0+45, 4); i0+=49;

    for (k=0;k<ngnss;k++) {
        if (i0+61>rtcm->havebit) return FALSE;
        cmi = getbitu(rtcm->buff, i0+60, 1); i0+=61;
        if (cmi) {
            svmask = (uint64_t)getbitu(rtcm->buff, i0, 8)<<32; i0+=8;
            svmask |= (uint64_t)getbitu(rtcm->buff, i0, 32); i0+=32;
            sigmask = getbitu(rtcm->buff, i0, 16); i0+=16;
            nsat = svmask2nsat(svmask);
            nsig = sigmask2nsig(sigmask);
            if (i0+nsat*nsig>rtcm->havebit) return FALSE;
            i0+=nsat*nsig;
        }
    }
    return TRUE;
}

static void output_cssr_oc(rtcm_t *rtcm, cssr_t *cssr, FILE *fp)
{
    int j, prn, sat[CSSR_MAX_SV], iode, nsat, gnss;
    ssr_t *ssr=NULL;

    if (fp == NULL) return;
    nsat = svmask2sat(cssr->svmask,sat);
    for (j=0;j<nsat;j++) {
        ssr = &rtcm->nav.ssr[sat[j]-1];
        gnss = sys2gnss(satsys(sat[j], &prn), NULL);
        iode = ssr->iode;
        if (j != 0) {
            fprintf(fp, ",,,,,, %d, %d, %d, ", gnss, prn, iode);
        } else {
            fprintf(fp, ", %d, %d, %d, ", gnss, prn, iode);
        }
        
        if (ssr->deph[0] != INVALID_VALUE) {
            fprintf(fp, "%f, ", (double)ssr->deph[0]);
        } else {
            fprintf(fp, "#N/A, ");
        }
        if (ssr->deph[1] != INVALID_VALUE) {
            fprintf(fp, "%f, ", (double)ssr->deph[1]);
        } else {
            fprintf(fp, "#N/A, ");
        }
        if (ssr->deph[2] != INVALID_VALUE) {
            fprintf(fp, "%f\n", (double)ssr->deph[2]);
        }else {
            fprintf(fp, "#N/A\n");
        }
    }
    if (nsat == 0) {
        fprintf(fp, "\n");
    }
}

/* decode orbit correction */
static int decode_cssr_oc(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, sys, iode;
    double udint;
    ssr_t *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    nsat = svmask2sat(cssr->svmask,sat);
    check_week_ref(rtcm, tow, ref_orbit);
    rtcm->time = gpst2time(rtcm->week_ref[ref_orbit], tow);

    trace(2,"decode_cssr_oc:   facility=%d tow=%d iod=%d\n", l6facility[chidx]+1, tow, iod);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (j = 0; j < MAXSAT; ++j) {
        ssr=&rtcm->nav.ssr[j];
        ssr->t0[0].sec = 0.0;
        ssr->t0[0].time = 0;
        ssr->udi[0]=0;
        ssr->iod[0]=0;
        ssr->update_oc=0;
        ssr->update=0;
        ssr->iode = 0;
        ssr->deph[0] = 0.0;
        ssr->deph[1] = 0.0;
        ssr->deph[2] = 0.0;
    }
    
    for (j=0;j<nsat;j++) {
        ssr = &rtcm->nav.ssr[sat[j]-1];
        sys = satsys(sat[j], NULL);

        if (sys == SYS_GAL) {
            iode = getbitu(rtcm->buff, i, 10); i+= 10; /* iode */
        } else {
            iode = getbitu(rtcm->buff, i, 8); i+= 8; /* iode */
        }

        /* delta radial/along-track/cross-track */
        ssr->deph[0] = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;
        ssr->deph[1] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
        ssr->deph[2] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
        if (ssr->deph[0] == INVALID_VALUE || ssr->deph[1] == INVALID_VALUE || ssr->deph[2] == INVALID_VALUE) {
            trace(3, "invalid orbit value: tow=%d, sat=%d, value=%d %d %d\n", 
                  tow, sat[j], ssr->deph[0], ssr->deph[1], ssr->deph[2]);
        }
        ssr->iode = iode;

        ssr->t0 [0]=rtcm->time;
        ssr->udi[0]=udint;
        ssr->iod[0]=cssr->iod;
        if (ssr->deph[0] == INVALID_VALUE || ssr->deph[1] == INVALID_VALUE || ssr->deph[2] == INVALID_VALUE) {
            ssr->deph[0] = INVALID_VALUE;
            ssr->deph[1] = INVALID_VALUE;
            ssr->deph[2] = INVALID_VALUE;
        }

        for (k=0;k<3;k++) ssr->ddeph[k]=0.0;

        ssr->update_oc=1;
        ssr->update=1;
        trace(4, "ssr orbit: prn=%2d, tow=%d, udi=%.1f, iod=%2d, orb=%f,%f,%f\n",sat[j],tow,
              udint,cssr->iod,ssr->deph[0],ssr->deph[1],ssr->deph[2]);
    }
    output_cssr_oc(rtcm, cssr, fp);

    check_cssr_changed_facility(cssr->l6facility);
    set_cssr_bank_orbit(rtcm->time, &rtcm->nav, 0);
    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is enough to decode the orbit correction message */
static int check_bit_width_oc(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int k,sat[CSSR_MAX_SV],nsat,prn;

    i0+=21;
    if (i0>rtcm->havebit) return FALSE;
    nsat = svmask2sat(cssr->svmask, sat);
    for (k=0;k<nsat;k++) {
        i0+=(satsys(sat[k],&prn)==SYS_GAL)?51:49;
        if (i0>rtcm->havebit) return FALSE;
    }
    return TRUE;
}

static void output_cssr_cc(rtcm_t *rtcm, cssr_t *cssr, FILE *fp)
{
    int j, tow, sat[CSSR_MAX_SV], nsat, prn, gnss, week;
    ssr_t *ssr=NULL;

    if (fp == NULL) return;
    tow = time2gpst(rtcm->time, &week);
    nsat = svmask2sat(cssr->svmask,sat);

    for (j=0;j<nsat;j++) {
        ssr = &rtcm->nav.ssr[sat[j]-1];
        gnss = sys2gnss(satsys(sat[j], &prn), NULL);
        if (j != 0) {
            fprintf(fp, ",,,,,, %d, %d, ", gnss, prn);
        } else {
            fprintf(fp, ", %d, %d, ", gnss, prn);
        }

        if (ssr->dclk[0] == INVALID_VALUE) {
            trace(3, "invalid clock value: tow=%d, sat=%d, value=%d\n", tow, sat[j], ssr->dclk[0]);
            fprintf(fp, "#N/A\n");
        } else {
            fprintf(fp, "%f\n", (double)ssr->dclk[0]);
        }
    }
    if (nsat == 0) {
        fprintf(fp, "\n");
    }
}

/* decode clock correction */
static int decode_cssr_cc(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, iod, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat;
    double udint;
    ssr_t *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_clock);
    rtcm->time = gpst2time(rtcm->week_ref[ref_clock], tow);
    nsat = svmask2sat(cssr->svmask,sat);

    trace(2,"decode_cssr_cc:   facility=%d tow=%d iod=%d\n", l6facility[chidx]+1, tow, iod);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (j = 0; j < MAXSAT; ++j) {
        ssr=&rtcm->nav.ssr[j];
        ssr->t0[1].sec = 0.0;
        ssr->t0[1].time = 0;
        ssr->udi[1]=0;
        ssr->iod[1]=0;
        ssr->update_cc=0;
        ssr->update=0;
        ssr->dclk[0] = 0.0;
    }
    
    for (j=0;j<nsat;j++) {
        ssr = &rtcm->nav.ssr[sat[j]-1];
        ssr->t0 [1]=rtcm->time;
        ssr->udi[1]=udint;
        ssr->iod[1]=cssr->iod;

        ssr->dclk[0] = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;
        if (ssr->dclk[0] == INVALID_VALUE) {
            trace(3, "invalid clock value: tow=%d, sat=%d, value=%d\n", tow, sat[j], ssr->dclk[0]);
        }
        ssr->dclk[1] = ssr->dclk[2] = 0.0;
        ssr->update_cc=1;
        ssr->update=1;
        trace(4, "ssr clock: prn=%2d, tow=%d, udi=%.1f, iod=%2d, clk=%f\n", sat[j], tow,
              udint, cssr->iod, ssr->dclk[0]);
    }
    output_cssr_cc(rtcm, cssr, fp);
    check_cssr_changed_facility(cssr->l6facility);
    set_cssr_bank_clock(rtcm->time, &rtcm->nav, 0);
    rtcm->nbit = i;
    return sync ? 0:10;
}

static int check_bit_width_cc(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int nsat;

    nsat = svmask2sat(cssr->svmask,NULL);
    return (i0+21+15*nsat<=rtcm->havebit);
}

static int sigmask2sig_p(int nsat, int *sat, uint16_t *sigmask, uint16_t *cellmask,
        int *nsig, int *sig)
{
    int j,k,id,sys,sys_p=-1,nsig_s=0,code[CSSR_MAX_SIG];

    for (j=0;j<nsat;j++) {
        sys = satsys(sat[j], NULL);
        if (sys != sys_p ){
            id = sys2gnss(sys, NULL);
            for (k=0,nsig_s=0;k<CSSR_MAX_SIG;k++) {
                if ((sigmask[id]>>(CSSR_MAX_SIG-1-k))&1) {
                    code[nsig_s] = k;
                    nsig_s++;
                }
            }
        }
        sys_p = sys;

        for (k=0, nsig[j]=0;k<nsig_s;k++) {
            if ((cellmask[j]>>(nsig_s-1-k))&1) {
                if (sig)
                    sig[j*CSSR_MAX_SIG+nsig[j]] = code[k];
                nsig[j]++;
            }
        }
    }

    return 1;
}


/* decode available signals from sigmask */
static int sigmask2sig(int nsat, int *sat, uint16_t *sigmask, uint16_t *cellmask,
        int *nsig, int *sig)
{
    int j,k,id,*codes=NULL,sys,sys_p=-1,ofst=0,nsig_s=0,code[CSSR_MAX_SIG],csize=0;
    const int codes_gps[]={
        CODE_L1C,CODE_L1P,CODE_L1W,CODE_L1S,CODE_L1L,CODE_L1X,
        CODE_L2S,CODE_L2L,CODE_L2X,CODE_L2P,CODE_L2W,
        CODE_L5I,CODE_L5Q,CODE_L5X
    };
    const int codes_glo[]={
        CODE_L1C,CODE_L1P,CODE_L2C,CODE_L2P,CODE_L3I,CODE_L3Q,CODE_L3X
    };
    const int codes_gal[]={
        CODE_L1B,CODE_L1C,CODE_L1X,CODE_L5I,CODE_L5Q,
        CODE_L5X,CODE_L7I,CODE_L7Q,CODE_L7X,CODE_L8I,CODE_L8Q,
        CODE_L8X
    };
    const int codes_qzs[]={
        CODE_L1C,CODE_L1S,CODE_L1L,CODE_L1X,CODE_L2S,CODE_L2L,CODE_L2X,
        CODE_L5I,CODE_L5Q,CODE_L5X
    };
    const int codes_bds[]={
        CODE_L2I,CODE_L2Q,CODE_L2X,
        CODE_L6I,CODE_L6Q,CODE_L6X,
        CODE_L7I,CODE_L7Q,CODE_L7X
    };
    const int codes_sbs[]={
        CODE_L1C,CODE_L5I,CODE_L5Q,CODE_L5X
    };

    for (j=0;j<nsat;j++,ofst+=nsig_s) {
        sys = satsys(sat[j], NULL);
        if (sys != sys_p) {
            id = sys2gnss(sys, NULL);
            ofst = 0;
            switch (sys) {
                case SYS_GPS: codes = (int *)codes_gps; csize=sizeof(codes_gps)/sizeof(int); break;
                case SYS_GLO: codes = (int *)codes_glo; csize=sizeof(codes_glo)/sizeof(int); break;
                case SYS_GAL: codes = (int *)codes_gal; csize=sizeof(codes_gal)/sizeof(int); break;
                case SYS_CMP: codes = (int *)codes_bds; csize=sizeof(codes_bds)/sizeof(int); break;
                case SYS_QZS: codes = (int *)codes_qzs; csize=sizeof(codes_qzs)/sizeof(int); break;
                case SYS_SBS: codes = (int *)codes_sbs; csize=sizeof(codes_sbs)/sizeof(int); break;
            }
            for (k=0,nsig_s=0;k<csize&&k<CSSR_MAX_SIG;k++) {
                if ((sigmask[id]>>(CSSR_MAX_SIG-1-k))&1) {
                    code[nsig_s] = codes[k];
                    nsig_s++;
                }
            }
        }
        sys_p = sys;

        for (k=0, nsig[j]=0;k<nsig_s;k++) {
            if ((cellmask[j]>>(nsig_s-1-k))&1) {
                if (sig)
                    sig[j*CSSR_MAX_SIG+nsig[j]] = code[k];
                nsig[j]++;
            }
        }
    }

    return 1;
}

static void output_cssr_cb(rtcm_t *rtcm, cssr_t *cssr, FILE *fp)   
{
    int k, j, sat[CSSR_MAX_SV], s, tow, lcnt, prn, nsat, gnss, week;
    ssr_t *ssr=NULL;
    int first = TRUE;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];

    if (fp == NULL) return;
    tow = time2gpst(rtcm->time, &week);
    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);

    for (k=lcnt=0;k<nsat;k++) {
        ssr = &rtcm->nav.ssr[sat[k]-1];
        for (j=0;j<nsig[k];j++) {
            gnss = sys2gnss(satsys(sat[k], &prn), NULL);
            if (first == FALSE) {
                fprintf(fp, ",,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
            } else {
                fprintf(fp, ", %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                first = FALSE;
            }
            s = sig[k*CSSR_MAX_SIG+j];
            /* code bias */
            if (ssr->cbias[s-1] == INVALID_VALUE) {
                trace(3, "invalid cb value: tow=%d, sat=%d, value=%.2f\n", tow, sat[k], ssr->cbias[s-1]);
                fprintf(fp, "#N/A\n");
            } else {
                fprintf(fp, "%f\n", ssr->cbias[s-1]);
            }
            ++lcnt;
            ssr->smode[j]=s;
        }
        ssr->nsig=nsig[k];
    }
    if (lcnt == 0) {
        fprintf(fp, "\n");
    }
}

/* decode code bias message */
static int decode_cssr_cb(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, s, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, lcnt;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];
    double udint;
    ssr_t *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_cbias);
    rtcm->time = gpst2time(rtcm->week_ref[ref_cbias], tow);
    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);

    trace(2,"decode_cssr_cb:   facility=%d tow=%d iod=%d\n", l6facility[chidx]+1, tow, iod);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (k = 0; k < MAXSAT; ++k) {
        ssr = &rtcm->nav.ssr[k];
        ssr->t0[4].sec = 0.0;
        ssr->t0[4].time = 0;
        ssr->udi[4]=0;
        ssr->iod[4]=0;
        ssr->update_cb=0;
        ssr->update=0;
        ssr->nsig=0;
        
        for (j = 0; j < MAXCODE; ++j) {
            /* code bias */
            ssr->cbias[j] = 0.0;
            ssr->smode[j] = 0;
        }
    }
    
    for (k=lcnt=0;k<nsat;k++) {
        ssr = &rtcm->nav.ssr[sat[k]-1];
        ssr->t0 [4]=rtcm->time;
        ssr->udi[4]=udint;
        ssr->iod[4]=cssr->iod;
        ssr->update_cb=1;
        ssr->update=1;

        for (j=0;j<nsig[k];j++) {
            s = sig[k*CSSR_MAX_SIG+j];
            /* code bias */
            ssr->cbias[s-1] = decode_sval(rtcm->buff, i, 11, 0.02); i+=11;
            if (ssr->cbias[s-1] == INVALID_VALUE) {
                trace(3, "invalid cb value: tow=%d, sat=%d, value=%d\n", tow, sat[k], ssr->cbias[s-1]);
            }
            trace(4, "ssr cbias: prn=%2d, tow=%d, udi=%.1f, iod=%2d, cbias=%f\n", sat[k], tow,
                  udint, cssr->iod, ssr->cbias[s-1]);
            ++lcnt;
            ssr->smode[j]=s;
        }
        ssr->nsig=nsig[k];
    }
    output_cssr_cb(rtcm, cssr, fp);
    check_cssr_changed_facility(cssr->l6facility);
    set_cssr_bank_cbias(rtcm->time, &rtcm->nav, 0, 0);
    rtcm->nbit = i;
    return sync ? 0:10;
}

static void output_cssr_pb(rtcm_t *rtcm, cssr_t *cssr, FILE *fp)
{
    int j, k, s, tow, sat[CSSR_MAX_SV], nsat, gnss, prn, week;
    int lcnt;
    int first = TRUE;
    ssr_t *ssr=NULL;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];

    if (fp == NULL) return;
    tow = time2gpst(rtcm->time, &week);
    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);

    for (k=lcnt=0;k<nsat;k++) {
        ssr = &rtcm->nav.ssr[sat[k]-1];
        gnss = sys2gnss(satsys(sat[k], &prn), NULL);

        for (j = 0; j < nsig[k]; j++) {
            s = sig[k*CSSR_MAX_SIG+j];
            if (first == FALSE) {
                fprintf(fp, ",,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
            } else {
                fprintf(fp, ", %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                first = FALSE;
            }
            if (ssr->pbias[s-1]  != INVALID_VALUE) {
                fprintf(fp, "%f, ", (double)ssr->pbias[s-1] );
                trace(4, "ssr pbias: prn=%2d, tow=%d, udi=%.1f, iod=%2d, pbias=%f\n", sat[k], tow,
                     ssr->udi[5], cssr->iod, ssr->pbias[s-1]);
            } else {
                fprintf(fp, "#N/A, ");
                trace(3, "invalid pb value: tow=%d, sat=%d, value=%d\n", tow, sat[k], ssr->pbias[s-1]);
            }
            fprintf(fp, "%d\n", ssr->discontinuity[s-1]);
            ++lcnt;
        }
    }
    if (lcnt == 0) {
        fprintf(fp, "\n");
    }

}

/* decode phase bias message */
static int decode_cssr_pb(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, s, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat;
    int lcnt;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];
    double udint;
    ssr_t *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_pbias);
    rtcm->time = gpst2time(rtcm->week_ref[ref_pbias], tow);
    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);

    trace(2,"decode_cssr_pb:   facility=%d tow=%d iod=%d\n", l6facility[chidx]+1, tow, iod);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (k = 0; k < MAXSAT; ++k) {
        ssr = &rtcm->nav.ssr[k];
        ssr->t0[5].sec = 0.0;
        ssr->t0[5].time = 0;
        ssr->udi[5]=0;
        ssr->iod[5]=0;
        ssr->update_pb=0;
        ssr->update=0;
        ssr->nsig=0;
        
        for (j = 0; j < MAXCODE; ++j) {
            /* phase bias */
            ssr->pbias[j] = 0.0;
            ssr->smode[j] = 0;
        }
    }


    for (k=lcnt=0;k<nsat;k++) {
        ssr = &rtcm->nav.ssr[sat[k]-1];
        ssr->t0 [5]=rtcm->time;
        ssr->udi[5]=udint;
        ssr->iod[5]=cssr->iod;
        ssr->update_pb=1;
        ssr->update=1;

        for (j=0;j<nsig[k];j++) {
            s = sig[k*CSSR_MAX_SIG+j];
            ssr->pbias[s-1] = decode_sval(rtcm->buff, i, 15,0.001); i+=15;
            ssr->discontinuity[s-1] = getbitu(rtcm->buff, i, 2); i+=2;
            ssr->smode[j]=s;
            ++lcnt;
        }
        ssr->nsig=nsig[k];
    }
    output_cssr_pb(rtcm, cssr, fp);
    check_cssr_changed_facility(cssr->l6facility);
    set_cssr_bank_pbias(rtcm->time, &rtcm->nav, 0, 0);

    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is sufficient to decode the code bias message */
static int check_bit_width_cb(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int nsig[CSSR_MAX_SV], nsig_total=0;
    int k, sat[CSSR_MAX_SV], nsat;

    if (rtcm->subtype!=CSSR_TYPE_CB && rtcm->subtype!=CSSR_TYPE_BIAS)
        return FALSE;
    
    nsat = svmask2sat(cssr->svmask, sat);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, NULL);

    for (k=0;k<nsat;k++) {
        nsig_total+=nsig[k];
    }
    return i0+21+nsig_total*11<=rtcm->havebit;
}

static int check_bit_width_pb(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int nsig[CSSR_MAX_SV], nsig_total=0;
    int k, sat[CSSR_MAX_SV], nsat;

    if (rtcm->subtype!=CSSR_TYPE_PB && rtcm->subtype!=CSSR_TYPE_BIAS)
        return FALSE;
    
    nsat = svmask2sat(cssr->svmask, sat);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, NULL);

    for (k=0;k<nsat;k++) {
        nsig_total+=nsig[k];
    }
    return i0+21+nsig_total*17<=rtcm->havebit;
}

static void output_cssr_bias(rtcm_t *rtcm, cssr_t *cssr, int cbflag, int pbflag, int netflag, int network, int netmask, FILE *fp)
{
    int j, k, s, sat[CSSR_MAX_SV], nsat;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];
    int lcnt, gnss, prn;
    ssr_t *ssr=NULL;

    if (fp == NULL) return;

    fprintf(fp, ", %d", cbflag);
    fprintf(fp, ", %d", pbflag);
    fprintf(fp, ", %d", netflag);
    if (netflag == 1) {
        fprintf(fp, ", %d", network);
    }

    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);
    if (netflag == 1) {
        fprintf(fp, ", 0x%08x", netmask);
    }

    for (k=lcnt=0;k<nsat;k++) {
        ssr = &rtcm->nav.ssr[sat[k]-1];
        if (!((netmask>>(nsat-1-k)) & 1)) {
            continue;
        }
        gnss = sys2gnss(satsys(sat[k], &prn), NULL);
        for (j = 0; j < nsig[k]; j++) {
            if (lcnt != 0) {
                if (netflag == 1) {
                    fprintf(fp, ",,,,,,,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                } else {
                    fprintf(fp, ",,,,,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                }
            } else {
                fprintf(fp, ", %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
            }
            s = sig[k*CSSR_MAX_SIG+j];
            if (cbflag == 1) { 
                /* code bias */
                if (ssr->cbias[s-1] != INVALID_VALUE) {
                    fprintf(fp, "%f, ", (double)ssr->cbias[s-1]);
                }
            } else {
                fprintf(fp, ", ");
            }
            if (pbflag == 1) {
                /* phase bias */
                if (ssr->pbias[s-1] != INVALID_VALUE) {
                    fprintf(fp, "%f, ", (double)ssr->pbias[s-1]);
                } else {
                    fprintf(fp, "#N/A, ");
                }
                fprintf(fp, "%d", ssr->discontinuity[s-1]);
            } else {
                fprintf(fp, ", ");
            }
            fprintf(fp, "\n");
            ++lcnt;
        }
    }
    if (lcnt == 0) {
        fprintf(fp, "\n");
    }

}

/* code bias correction */
static int decode_cssr_bias(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, s, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];
    int cbflag, pbflag, netflag, network, netmask, lcnt;
    double udint;
    ssr_t *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_bias);
    rtcm->time = gpst2time(rtcm->week_ref[ref_bias], tow);

    cbflag = getbitu(rtcm->buff, i, 1); i += 1;
    pbflag = getbitu(rtcm->buff, i, 1); i += 1;
    netflag = getbitu(rtcm->buff, i, 1); i += 1;
    network = getbitu(rtcm->buff, i, (netflag ? 5: 0)); i += (netflag ? 5: 0);

    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);
    netmask = getbitu(rtcm->buff, i, (netflag ? nsat: 0)); i += (netflag ? nsat: 0);

    trace(2,"decode_cssr_bias: facility=%d tow=%d iod=%d net=%2d mask=0x%x flag=%d %d %d\n", l6facility[chidx]+1, tow, iod, network, netmask, cbflag, pbflag, netflag);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }
    
    for (k = 0; k < MAXSAT; ++k) {
        ssr = &rtcm->nav.ssr[k];
        
        if (rtcm->subtype == CSSR_TYPE_BIAS && cbflag == 1) {
            ssr->t0[4].sec = 0.0;
            ssr->t0[4].time = 0;
            ssr->udi[4]=0;
            ssr->iod[4]=0;
            ssr->update_cb=0;
            ssr->update=0;
            ssr->nsig=0;
        }
        if (rtcm->subtype == CSSR_TYPE_BIAS && pbflag == 1) {
            ssr->t0[5].sec = 0.0;
            ssr->t0[5].time = 0;
            ssr->udi[5]=0;
            ssr->iod[5]=0;
            ssr->update_pb=0;
            ssr->update=0;
            ssr->nsig=0;
        }
        
        for (j = 0; j < MAXCODE; ++j) {
            if (rtcm->subtype == CSSR_TYPE_BIAS && cbflag == 1) {
                /* code bias */
                ssr->smode[j] = 0;
                ssr->cbias[j] = 0.0;
            }
            if (rtcm->subtype == CSSR_TYPE_BIAS && pbflag == 1) {
                /* phase bias */
                ssr->smode[j] = 0;
                ssr->pbias[j] = 0.0;
                ssr->discontinuity[j] = 0;
            }
        }
    }
    
    for (k=lcnt=0;k<nsat;k++) {
        if (netflag && !((netmask>>(nsat-1-k)) & 1)) {
            continue;
        }
        ssr = &rtcm->nav.ssr[sat[k]-1];

        if (rtcm->subtype == CSSR_TYPE_BIAS && cbflag == 1) {
            ssr->t0 [4]=rtcm->time;
            ssr->udi[4]=udint;
            ssr->iod[4]=cssr->iod;
            ssr->update_cb=1;
            ssr->update=1;
        }
        if (rtcm->subtype == CSSR_TYPE_BIAS && pbflag == 1) {
            ssr->t0 [5]=rtcm->time;
            ssr->udi[5]=udint;
            ssr->iod[5]=cssr->iod;
            ssr->update_pb=1;
            ssr->update=1;
        }

        for (j=0;j<nsig[k];j++) {
            s = sig[k*CSSR_MAX_SIG+j];
            if (cbflag == 1) { /* code bias */
                ssr->cbias[s-1] = decode_sval(rtcm->buff, i, 11, 0.02); i+=11;
                if (ssr->cbias[s-1] == INVALID_VALUE) {
                    trace(3, "invalid cb value: tow=%d, sat=%d, value=%d\n", tow, sat[k], ssr->cbias[s-1]);
                }
                trace(4, "ssr cbias: network=%d, prn=%2d, tow=%d, udi=%.1f, iod=%2d, s=%d, cbias=%f\n",
                    network, sat[k], tow, udint, cssr->iod, s, ssr->cbias[s-1]);
            }
            if (pbflag == 1) { /* phase bias */
                /* phase bias */
                ssr->pbias[s-1] = decode_sval(rtcm->buff, i, 15, 0.001); i+=15;
                if (ssr->pbias[s-1] == INVALID_VALUE) {
                    trace(3, "invalid pb value: tow=%d, sat=%d, value=%d\n", tow, sat[k], ssr->pbias[s-1]);
                }
                ssr->discontinuity[s-1] = getbitu(rtcm->buff, i, 2); i+= 2;
                trace(4, "ssr pbias: network=%d, prn=%2d, tow=%d, udi=%.1f, iod=%2d, s=%d, pbias=%f\n",
                    network, sat[k], tow, udint, cssr->iod, s, ssr->pbias[s-1]);
            }
            ++lcnt;
            ssr->smode[j]=s;
        }
        ssr->nsig=nsig[k];
    }
    output_cssr_bias(rtcm, cssr, cbflag, pbflag, netflag, network, netmask, fp);
    if (cbflag == 1) {
        check_cssr_changed_facility(cssr->l6facility);
        set_cssr_bank_cbias(rtcm->time, &rtcm->nav, (netflag ? network: 0), cssr->iod);
    }
    if (pbflag == 1) {
        check_cssr_changed_facility(cssr->l6facility);
        set_cssr_bank_pbias(rtcm->time, &rtcm->nav, (netflag ? network: 0), cssr->iod);
    }
    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is sufficient to decode the bias message */
static int check_bit_width_bias(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int j,k,nsat,slen=0,cbflag,pbflag,netflag,netmask=0;
    int sat[CSSR_MAX_SV], nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG];

    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);
    if (i0+24>rtcm->havebit) return FALSE;

    i0+=21;
    cbflag = getbitu(rtcm->buff, i0, 1); i0+=1;
    pbflag = getbitu(rtcm->buff, i0, 1); i0+=1;
    netflag = getbitu(rtcm->buff, i0, 1); i0+=1;

    if (netflag) {
        if (i0+5+nsat>rtcm->havebit) return FALSE;
        i0+=5;
        netmask = getbitu(rtcm->buff, i0, nsat); i0+=nsat;
    }

    if (cbflag) slen+=11;
    if (pbflag) slen+=17;

    for (k=0;k<nsat;k++) {
        if (netflag && !((netmask>>(nsat-1-k))&1)) continue;
        for (j = 0; j < nsig[k]; j++) {
            if (i0+slen>rtcm->havebit) return FALSE;
            i0 += slen;
        }
    }
    return TRUE;
}

static void output_cssr_ura(rtcm_t *rtcm, cssr_t *cssr, FILE *fp)
{
    int j, sat[CSSR_MAX_SV], nsat, gnss, prn;
    ssr_t *ssr = NULL;

    if (fp == NULL) return;
    nsat = svmask2sat(cssr->svmask,sat);
    for (j=0;j<nsat;j++) {
        ssr = &rtcm->nav.ssr[sat[j]-1];
        gnss = sys2gnss(satsys(sat[j], &prn), NULL);
        if (j != 0) {
            fprintf(fp, ",,,,,, %d, %d, 0x%02x\n", gnss, prn, ssr->ura);
        } else {
            fprintf(fp, ", %d, %d, 0x%02x\n", gnss, prn, ssr->ura);
        }
    }
    if (nsat == 0) {
        fprintf(fp, "\n");
    }
}

/* decode ura correction */
static int decode_cssr_ura(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, iod, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat;
    double udint;
    ssr_t *ssr = NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_ura);
    rtcm->time = gpst2time(rtcm->week_ref[ref_ura], tow);
    nsat = svmask2sat(cssr->svmask,sat);

    trace(3,"decode_cssr_ura:  facility=%d tow=%d iod=%d\n", l6facility[chidx]+1, tow, iod);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (j=0;j<nsat;j++) {
        ssr = &rtcm->nav.ssr[sat[j]-1];
        ssr->t0 [3]=rtcm->time;
        ssr->udi[3]=udint;
        ssr->iod[3]=iod;
        ssr->ura = getbitu(rtcm->buff, i, 6); i+= 6; /* ssr ura */
        ssr->update_ura=1;
        ssr->update=1;
    }
    output_cssr_ura(rtcm, cssr, fp);

    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is sufficient to decode the ura message */
static int check_bit_width_ura(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int nsat;

    nsat = svmask2sat(cssr->svmask, NULL);
    return i0+21+6*nsat<=rtcm->havebit;
}

static void output_cssr_stec(rtcm_t *rtcm, cssr_t *cssr, int *qual_stec, int inet, FILE *fp)
{
    int j, s, sat[CSSR_MAX_SV], nsat, gnss, prn;
    ssrion_t *ssr_ion;

    if (fp == NULL) return;
    ssr_ion = &rtcm->ssr_ion[inet];

    nsat = svmask2sat(cssr->svmask, sat);

    fprintf(fp, ", %d, ", cssr->opt.stec_type);
    fprintf(fp, "%d, ", inet);

    fprintf(fp, "0x%02x", (uint32_t)(cssr->net_svmask[inet]>>32));
    fprintf(fp, "%08lx", cssr->net_svmask[inet]&0xffffffff);

    for (j=0,s=0;j<nsat;j++) {
        if ((cssr->net_svmask[inet]>>(nsat-1-j))&1) {
            gnss = sys2gnss(satsys(sat[j], &prn), NULL);
            if (s != 0) {
                fprintf(fp, ",,,,,,,,, %d, %d, ", gnss, prn);
            } else {
                fprintf(fp, ", %d, %d, ", gnss, prn);
            }
            fprintf(fp, "0x%02x, ", qual_stec[j]);
            
            if (ssr_ion->stec.a[s][0] != INVALID_VALUE) {
                fprintf(fp, "%f", (double)ssr_ion->stec.a[s][0]);
            } else {
                fprintf(fp, "#N/A");
            }
            if (cssr->opt.stec_type > 0) {
                if (ssr_ion->stec.a[s][1] != INVALID_VALUE) {
                    fprintf(fp, ", %f, ", (double)ssr_ion->stec.a[s][1]);
                } else {
                    fprintf(fp, ", #N/A, ");
                }

                if (ssr_ion->stec.a[s][2] != INVALID_VALUE) {
                    fprintf(fp, "%f", (double)ssr_ion->stec.a[s][2]);
                } else {
                    fprintf(fp, "#N/A");
                }
            }
            if (cssr->opt.stec_type > 1) {
                if (ssr_ion->stec.a[s][3] != INVALID_VALUE) {
                    fprintf(fp, ", %f", (double)ssr_ion->stec.a[s][3]);
                } else {
                    fprintf(fp, ", #N/A");
                }
            }
            fprintf(fp, "\n");
            s++;
        }
    }
    if (s == 0) {
        fprintf(fp, "\n");
    }
}

/* decode stec correction */
static int decode_cssr_stec(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, s, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, inet, a, b, qual_stec[CSSR_MAX_SV];
    double udint;
    ssrgp_t *ssrg;
    ssrion_t *ssr_ion;
    
    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_stec);
    rtcm->time = gpst2time(rtcm->week_ref[ref_stec], tow);
    nsat = svmask2sat(cssr->svmask, sat);
    cssr->opt.stec_type = getbitu(rtcm->buff, i, 2); i+= 2; /* stec correction type */
    inet = getbitu(rtcm->buff, i, 5); i+= 5; /* network id */

    trace(2,"decode_cssr_stec: facility=%d tow=%d iod=%d net=%2d type=%d\n", l6facility[chidx]+1, tow, iod, inet, cssr->opt.stec_type);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    ssrg = &rtcm->ssrg[inet];
    ssr_ion = &rtcm->ssr_ion[inet];
    
    cssr->net_svmask[inet] = getbitu(rtcm->buff, i, nsat); i+= nsat; /* stec correction type */
    trace(4, "decode_cssr_stec: mask=0x%x\n", cssr->net_svmask[inet]);

    ssrg->t0 = rtcm->time;
    ssrg->udi = udint;
    ssrg->iod = iod;

    for (j=0,s=0;j<nsat;j++) {
        if ((cssr->net_svmask[inet]>>(nsat-1-j))&1) {
            ssr_ion->stec.sat[s] = sat[j];
            a = getbitu(rtcm->buff, i, 3); i+= 3;
            b = getbitu(rtcm->buff, i, 3); i+= 3;
            ssr_ion->stec.quality[s] = decode_cssr_quality_stec(a,b);
            qual_stec[j] = (a<<3)|b;
            for (k=0;k<4;k++) ssr_ion->stec.a[s][k] = 0.0;

            ssr_ion->stec.a[s][0] = decode_sval(rtcm->buff, i, 14, 0.05); i+=14;
            if (cssr->opt.stec_type > 0) {
                ssr_ion->stec.a[s][1] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
                ssr_ion->stec.a[s][2] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
            }
            if (cssr->opt.stec_type > 1) {
                ssr_ion->stec.a[s][3] = decode_sval(rtcm->buff, i, 10, 0.02); i+=10;
            }
            trace(4, "decode_cssr_stec: tow=%d, sat=%d\n", tow, sat[j]);
            s++;
        }
    }
    output_cssr_stec(rtcm, cssr, qual_stec, inet, fp);
    ssr_ion->stec.network = inet;
    ssr_ion->stec.nsat = s;
    ssrg->update = 1;

    rtcm->nbit = i;
    trace(3, "decode_cssr_stec(): tow=%d, net=%d, bits=%d\n", tow, inet, i - i0);
    return sync ? 0:10;
}

/* check if the buffer length is sufficient to decode the stec message */
static int check_bit_width_stec(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int j,sat[CSSR_MAX_SV],nsat,stec_type,slen=0,nsat_local=0;
    uint64_t net_svmask;
    const int slen_t[4] = {20,44,54,70};

    nsat = svmask2sat(cssr->svmask, sat);
    if (i0+28+nsat>rtcm->havebit) return FALSE;

    i0+=21;
    stec_type = getbitu(rtcm->buff, i0, 2); i0+=2;
    i0+=5;
    net_svmask = getbitu(rtcm->buff, i0, nsat); i0+=nsat;

    slen = slen_t[stec_type];

    for (j=0;j<nsat;j++) { /* number of local satellites */
        if ((net_svmask>>(nsat-1-j))&1) {
            nsat_local++;
        }
    }
    return i0+nsat_local*slen<=rtcm->havebit;
}

static void output_cssr_grid(rtcm_t *rtcm, cssr_t *cssr, int trop_type, int sz_idx, int trop_qual, int inet, double dstec_[][MAXSAT], FILE *fp)
{
    int j, k, ii, s, sat[CSSR_MAX_SV], nsat;
    int gnss, prn;
    ssrgp_t *ssrg;

    if (fp == NULL) return;
    nsat = svmask2sat(cssr->svmask,sat);

    fprintf(fp, ", %d, ", trop_type);
    fprintf(fp, "%d, ", sz_idx);
    fprintf(fp, "%d, ", inet);

    ssrg = &rtcm->ssrg[inet];

    fprintf(fp, "0x%02x", (uint32_t)(cssr->net_svmask[inet]>>32));
    fprintf(fp, "%08lx, ", cssr->net_svmask[inet]&0xffffffff);
    fprintf(fp, "0x%02x, ", trop_qual);
    fprintf(fp, "%d", ssrg->ngp);

    for (j=0;j<ssrg->ngp;j++) {
        if (j != 0) {
            fprintf(fp, ",,,,,,,,,,,, %d, ", j+1);
        } else {
            fprintf(fp, ", %d, ", j+1);
        }
        switch (trop_type) {
            case 0: break;
            case 1:
                if (ssrg->trop_wet[j] != INVALID_VALUE || ssrg->trop_total[j] != INVALID_VALUE) {
                    fprintf(fp, "%f, ", (double)ssrg->trop_total[j]-(double)ssrg->trop_wet[j]);
                    fprintf(fp, "%f", (double)ssrg->trop_wet[j]);
                } else {
                    fprintf(fp, "#N/A, ");
                    fprintf(fp, "#N/A");
                }
                break;
        }

        for (k=0,s=0,ii=0;k<nsat;k++) {
            if ((cssr->net_svmask[inet]>>(nsat-1-k)) & 1) {
                gnss = sys2gnss(satsys(sat[k], &prn), NULL);
                if (ii != 0) {
                    fprintf(fp, ",,,,,,,,,,,,,,, %d, %d, ", gnss, prn);
                } else {
                    fprintf(fp, ", %d, %d, ", gnss, prn);
                }
                if (ssrg->stec[j][s] == INVALID_VALUE) {
                    fprintf(fp, "#N/A\n");
                } else {
                    fprintf(fp, "%f\n", (double)dstec_[j][s]);
                }
                ++ii;
                ++s;
            }
        }
        if (ii == 0) {
            fprintf(fp, "\n");
        }
    }
    if (ssrg->ngp == 0) {
        fprintf(fp, "\n");
    }
}

/* decode grid correction */
static int decode_cssr_grid(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, ii, s, sync, iod, tow, ngnss, sat[CSSR_MAX_SV], trop_qual, nsat, sz;
    int trop_type, sz_idx, inet, a,b, hs, wet;
    double udint, stec0, dlat, dlon, dstec, dstec_[RTCM_SSR_MAX_GP][MAXSAT];
    nav_t *nav = &rtcm->nav;
    ssrgp_t *ssrg;
    ssrion_t *ssr_ion;
    int valid;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_grid);
    rtcm->time = gpst2time(rtcm->week_ref[ref_grid], tow);
    nsat = svmask2sat(cssr->svmask,sat);

    trop_type = getbitu(rtcm->buff, i, 2); i+=2;  /* troposphere correction type */
    sz_idx = getbitu(rtcm->buff, i, 1); i++; /* stec range */
    inet = getbitu(rtcm->buff, i, 5); i+=5; /* network id */

    ssrg = &rtcm->ssrg[inet];
    ssr_ion = &rtcm->ssr_ion[inet];

    cssr->net_svmask[inet] = getbitu(rtcm->buff, i, nsat); i+= nsat; /* stec correction type */
    a = getbitu(rtcm->buff, i, 3); i+= 3;
    b = getbitu(rtcm->buff, i, 3); i+= 3;
    trop_qual = (a<<3)|b;
    ssrg->ngp = getbitu(rtcm->buff, i, 6); i+= 6;
    ssrg->t0 = rtcm->time;
    ssrg->quality = decode_cssr_quality_trop(a,b);
    ssrg->network = inet;

    trace(2,"decode_cssr_grid: facility=%d tow=%d iod=%d net=%2d trop=%d sz_idx=%d ngp=%d\n", l6facility[chidx]+1, tow, iod, inet, trop_type, sz_idx, ssrg->ngp);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (j=0;j<RTCM_SSR_MAX_GP;j++) {
        ssrg->gp[j].pos[0] = 0.0;
        ssrg->gp[j].pos[1] = 0.0;
        ssrg->gp[j].pos[2] = 0.0;
        ssrg->gp[j].network = 0;
        ssrg->gp[j].update = 0;
        ssrg->nsv[j] = 0;
    }
    
    for (j=0;j<ssrg->ngp;j++) {
        ssrg->gp[j].pos[0] = clas_grid[inet][j][0]*D2R;
        ssrg->gp[j].pos[1] = clas_grid[inet][j][1]*D2R;
        ssrg->gp[j].pos[2] = clas_grid[inet][j][2];
        ssrg->gp[j].network = inet;
        ssrg->gp[j].update = 1;

        trace(4,"gp check:pos=%f,%f,%f,%d,%d\n",ssrg->gp[j].pos[0]*R2D,
              ssrg->gp[j].pos[1]*R2D,ssrg->gp[j].pos[2],inet,ssrg->ngp);
    }
    sz = (sz_idx)?16:7;

    for (j=0;j<ssrg->ngp;j++) {
        valid=1;
        switch (trop_type) {
            case 0: break;
            case 1:
                hs = getbits(rtcm->buff,i,9); i+=9;
                wet = getbits(rtcm->buff,i,8); i+= 8;
                if (hs==(-P2_S9_MAX-1)) {
                    trace(2, "trop(hs) is invalid: tow=%d, inet=%d, grid=%d, hs=%d\n",
                        tow, inet, j, hs);
                    valid=0;
                }
                if (wet==(-P2_S8_MAX-1)) {
                    trace(2, "trop(wet) is invalid: tow=%d, inet=%d, grid=%d, wet=%d\n",
                        tow, inet, j, wet);
                    valid=0;
                }
                if (valid == 1) {
                    ssrg->trop_wet[j] = wet*0.004+0.252;
                    ssrg->trop_total[j] = (hs+wet)*0.004+0.252+CSSR_TROP_HS_REF;
                } else {
                    ssrg->trop_wet[j] = INVALID_VALUE;
                    ssrg->trop_total[j] = INVALID_VALUE;
                }
                trace(4, "decode_cssr_grid: grid=%d, total=%.3f, wet=%.3f\n", j, ssrg->trop_total[j], ssrg->trop_wet[j]);
                break;
        }

        dlat = (ssrg->gp[j].pos[0] - ssrg->gp[0].pos[0])*R2D;
        dlon = (ssrg->gp[j].pos[1] - ssrg->gp[0].pos[1])*R2D;

        for (k=0,s=0,ii=0;k<nsat;k++) {
            if ((cssr->net_svmask[inet]>>(nsat-1-k)) & 1) {
                dstec_[j][s] = dstec = decode_sval(rtcm->buff, i, sz, 0.04); i+=sz;
                stec0 = ssr_ion->stec.a[ii][0] + ssr_ion->stec.a[ii][1]*dlat +
                        ssr_ion->stec.a[ii][2]*dlon +
                        ssr_ion->stec.a[ii][3]*dlat*dlon;
                if (dstec == INVALID_VALUE) {
                    trace(2, "dstec is invalid: tow=%d, inet=%d, grid=%d, sat=%d, dstec=%f\n",
                        tow, inet, j, sat[k], dstec);
                    ssrg->stec[j][s] = INVALID_VALUE;
                } else {
                    ssrg->stec[j][s] = dstec;
                    ssrg->stec[j][s] += stec0;
                }
                ssrg->stec0[j][s] = stec0;
                ssrg->sat[j][s] = sat[k];
                ++ii;
                ++s;
            }
        }
        ssrg->nsv[j] = s;
    }
    output_cssr_grid(rtcm, cssr, trop_type, sz_idx, trop_qual, inet, dstec_, fp);
    ssrg->update = 1;
    nav->updateac=1;

    check_cssr_changed_facility(cssr->l6facility);
    set_cssr_latest_trop(ssrg->t0, ssrg, inet);
    set_cssr_bank_trop(ssrg->t0, ssrg, inet);
    
    rtcm->nbit = i;
    trace(3, "decode_cssr_grid(): tow=%d, net=%d, bits=%d\n", tow, inet, i - i0);
    return sync ? 0:10;
}

/* check if the buffer length is sufficient to decode the grid message */
static int check_bit_width_grid(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int k,nsat,trop_type,ngp,sz_trop,sz_idx,sz_stec,nsat_local=0;
    uint64_t net_svmask;

    nsat = svmask2sat(cssr->svmask, NULL);
    if (i0+41+nsat>rtcm->havebit) return FALSE;
    i0 += 21;
    trop_type = getbitu(rtcm->buff, i0, 2); i0+=2;
    sz_idx = getbitu(rtcm->buff, i0, 1); i0++;
    i0+=5; /* network id */
    net_svmask = getbitu(rtcm->buff, i0, nsat); i0+=nsat;
    i0+=6; /* trop quality indicator */
    ngp = getbitu(rtcm->buff, i0, 6); i0+=6;

    sz_trop = (trop_type==0) ? 0:17;
    sz_stec = (sz_idx==0) ? 7:16;

    for (k=0;k<nsat;k++) {
        if ((net_svmask>>(nsat-1-k))&1) {
            nsat_local++;
        }
    }

    return i0+ngp*(sz_trop+nsat_local*sz_stec)<=rtcm->havebit;
}

static void output_cssr_combo(rtcm_t *rtcm, cssr_t *cssr, int flg_orbit, int flg_clock, int flg_net, int netid, int net_svmask, FILE *fp)
{
    int j, sat[CSSR_MAX_SV], nsat;
    int s, gnss, prn;
    ssr_t *ssr=NULL;

    if (fp == NULL) return;

    nsat = svmask2sat(cssr->svmask, sat);

    fprintf(fp, ", %d", flg_orbit);
    fprintf(fp, ", %d", flg_clock);
    fprintf(fp, ", %d", flg_net);
    fprintf(fp, ", %d", netid);
    fprintf(fp, ", 0x%04x", net_svmask & 0xffff);

    for (j = s = 0; j < nsat; ++j) {
        if ((net_svmask >> (nsat - 1 - j)) & 1) {
            ssr=&rtcm->nav.ssr[sat[j]-1];
            gnss = sys2gnss(satsys(sat[j], &prn), NULL);
            if (s != 0) {
                fprintf(fp, ",,,,,,,,,,, %d, %d", gnss, prn);
            } else {
                fprintf(fp, ", %d, %d", gnss, prn);
            }

            if (flg_orbit == 1) {
                fprintf(fp, ", %d", ssr->iode);
                /* delta radial,along-track,cross-track */
                if (ssr->deph[0] != INVALID_VALUE) {
                    fprintf(fp, ", %f", (double)ssr->deph[0]);
                } else {
                    fprintf(fp, ", #N/A");
                }
                if (ssr->deph[1] != INVALID_VALUE) {
                    fprintf(fp, ", %f", (double)ssr->deph[1]);
                } else {
                    fprintf(fp, ", #N/A");
                }
                if (ssr->deph[2] != INVALID_VALUE) {
                    fprintf(fp, ", %f", (double)ssr->deph[2]);
                } else {
                    fprintf(fp, ", #N/A");
                }
            } else {
                fprintf(fp, ",,,,");
            }
            if (flg_clock == 1) {
                if (ssr->dclk[0] == INVALID_VALUE) {
                    fprintf(fp, ", #N/A");
                } else {
                    fprintf(fp, ", %f", (double)ssr->dclk[0]);
                }
            }
            fprintf(fp, "\n");
            ++s;
        }
    }
    if (s == 0) {
        fprintf(fp, "\n");
    }
}

/* decode orbit/clock combination message */
static int decode_cssr_combo(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, sync, iod, tow, ngnss, sat[CSSR_MAX_SV], nsat, iode;
    int flg_orbit, flg_clock, flg_net, net_svmask, netid, s;
    double udint;
    ssr_t *ssr=NULL;
    
    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_combined);
    rtcm->time = gpst2time(rtcm->week_ref[ref_combined], tow);
    nsat = svmask2sat(cssr->svmask, sat);
    
    flg_orbit = getbitu(rtcm->buff, i, 1); i+=1;
    flg_clock = getbitu(rtcm->buff, i, 1); i+=1;
    flg_net = getbitu(rtcm->buff, i, 1); i+=1;
    netid = getbitu(rtcm->buff, i, (flg_net ? 5: 0)); i+=(flg_net ? 5: 0);
    net_svmask = getbitu(rtcm->buff, i, (flg_net ? nsat: 0)); i+=(flg_net ? nsat: 0);
    
    trace(2, "decode_cssr_combo:facility=%d tow=%d iod=%d net=%2d flag=%d %d %d\n", l6facility[chidx]+1, tow, iod, netid, flg_orbit, flg_clock, flg_net);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }
    
    for (j = 0; j < MAXSAT; ++j) {
        ssr=&rtcm->nav.ssr[j];
        if (flg_orbit == 1) {
            ssr->t0[0].sec = 0.0;
            ssr->t0[0].time = 0;
            ssr->udi[0]=0;
            ssr->iod[0]=0;
            ssr->update_oc = 0;
            ssr->update=0;
            ssr->iode = 0;
            ssr->deph[0] = 0.0;
            ssr->deph[1] = 0.0;
            ssr->deph[2] = 0.0;
        }
        if (flg_clock == 1) {
            ssr->t0[1].sec = 0.0;
            ssr->t0[1].time = 0;
            ssr->udi[1]=0;
            ssr->iod[1]=0;
            ssr->update_cc = 0;
            ssr->update=0;
            ssr->dclk[0] = 0.0;
        }
    }
    
    for (j = s = 0; j < nsat; ++j) {
        if (flg_net && !((net_svmask >> (nsat - 1 - j)) & 1)) {
            continue;
        }
        ssr=&rtcm->nav.ssr[sat[j]-1];
        
        if (flg_orbit == 1) {
            if (satsys(sat[j], NULL) == SYS_GAL) {
                iode = getbitu(rtcm->buff, i, 10); i += 10; /* iode */
            } else {
                iode = getbitu(rtcm->buff, i,  8); i +=  8; /* iode */
            }
            
            /* delta radial,along-track,cross-track */
            ssr->deph[0] = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;
            ssr->deph[1] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
            ssr->deph[2] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
            if (ssr->deph[0] == INVALID_VALUE || ssr->deph[1] == INVALID_VALUE || ssr->deph[2] == INVALID_VALUE) {
                trace(3, "invalid orbit value: tow=%d, sat=%d, value=%d %d %d\n", 
                      tow, sat[j], ssr->deph[0], ssr->deph[1], ssr->deph[2]);
            }
            ssr->iode = iode;
            
            ssr->t0 [0] = rtcm->time;
            ssr->udi[0] = udint;
            ssr->iod[0] = cssr->iod;
            
            for (k=0;k<3;k++) ssr->ddeph[k]=0.0;
            
            ssr->update_oc = 1;
            ssr->update=1;
            
            trace(4, "combined orbit: network=%d, tow=%d, sat=%d, iode=%d, deph=%f, %f, %f\n", netid, tow,
                sat[j], ssr->iode, ssr->deph[0], ssr->deph[1], ssr->deph[2]);
            if (ssr->deph[0] == INVALID_VALUE || ssr->deph[1] == INVALID_VALUE || ssr->deph[2] == INVALID_VALUE) {
                trace(3, "invalid orbit value: tow=%d, sat=%d, value=%f %f %f\n", tow, sat[j], ssr->deph[0], ssr->deph[1], ssr->deph[2]);
                ssr->deph[0] = INVALID_VALUE;
                ssr->deph[1] = INVALID_VALUE;
                ssr->deph[2] = INVALID_VALUE;
            }
        }
        if (flg_clock == 1) {
            ssr->dclk[0] = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;
            if (ssr->dclk[0] == INVALID_VALUE) {
                trace(3, "invalid clock value: tow=%d, sat=%d, value=%d\n", tow, sat[j], ssr->dclk[0]);
            }
            ssr->dclk[1] = ssr->dclk[2] = 0.0;
            ssr->t0 [1] = rtcm->time;
            ssr->udi[1] = udint;
            ssr->iod[1] = cssr->iod;
            ssr->update_cc = 1;
            ssr->update = 1;
            trace(4, "combined clock: network=%d, tow=%d, sat=%d, dclk=%f\n", netid, tow,
                sat[j], ssr->dclk[0]);
        }
        ++s;
    }
    output_cssr_combo(rtcm, cssr, flg_orbit, flg_clock, flg_net, netid, net_svmask, fp);
    check_cssr_changed_facility(cssr->l6facility);
    if (flg_net == 1) {
        cssrObject[chidx].separation |= (1 << (netid - 1));
    }

    if (flg_orbit == 1) {
        set_cssr_bank_orbit(rtcm->time, &rtcm->nav, (flg_net ? netid: 0));
    }
    if (flg_clock == 1) {
        set_cssr_bank_clock(rtcm->time, &rtcm->nav, (flg_net ? netid: 0));
    }

    rtcm->nbit = i;
    return sync ? 0: 10;
}

static int check_bit_width_combo(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int sat[CSSR_MAX_SV],nsat,j,flg_orbit,flg_clock,flg_net,sz;
    uint64_t net_svmask=0;

    nsat = svmask2sat(cssr->svmask,sat);
    if (i0+24>rtcm->havebit) return FALSE;

    i0+=21;
    flg_orbit = getbitu(rtcm->buff, i0, 1); i0+=1;
    flg_clock = getbitu(rtcm->buff, i0, 1); i0+=1;
    flg_net = getbitu(rtcm->buff, i0, 1); i0+=1;

    if (flg_net) {
        if (i0+5+nsat>rtcm->havebit) return FALSE;
        i0+=5; /* network id */
        net_svmask = getbitu(rtcm->buff, i0, nsat); i0+=nsat;
    }
    
    for (j=0;j<nsat;j++) {
        if (flg_net && !((net_svmask>>(nsat-1-j))&1)) {
            continue;
        }
        if (flg_orbit) {
            sz = (satsys(sat[j],NULL)==SYS_GAL) ? 10:8;
            i0+=sz+41;
            if (i0>rtcm->havebit) return FALSE;
        }
        if (flg_clock) {
            if ((i0+=15)>rtcm->havebit) return FALSE;
        }
    }
    return TRUE;
}

static void output_cssr_atmos(rtcm_t *rtcm, cssr_t *cssr, int trop_ctype, int stec_ctype, int inet, int trop_qual, int trop_type,
                             int sz_idx_t, int *sz_idx_s, double trop_ofst, double *ct, double ci_[][6], int *stec_qual, int *stec_type_, double *total, 
                             double *wet, double stec[][CSSR_MAX_SV], FILE *fp1, FILE *fp2)
{
    int j, k, s, gnss, sat[CSSR_MAX_SV], nsat, prn;
    const int ct_num[3]={1, 3, 4}, ci_num[4]={1, 3, 4, 6};
    ssrgp_t *ssrg;


    if (fp1 == NULL || fp2 == NULL) return;
    ssrg = &rtcm->ssrg[inet];
    nsat = svmask2sat(cssr->svmask, sat);

    fprintf(fp1, ", %d", trop_ctype);
    fprintf(fp1, ", %d", stec_ctype);
    fprintf(fp1, ", %d", inet);
    fprintf(fp2, ", %d", inet);
    fprintf(fp1, ", %d", ssrg->ngp);
    if (trop_ctype != 0) {
        fprintf(fp1, ", 0x%02x", trop_qual);
    } else {
        fprintf(fp1, ","); 
    }

    if ((trop_ctype&0x01) == 0x01) {
        fprintf(fp1, ", %d", trop_type);
        for (k=0; k<4; k++) {
            if (k < ct_num[trop_type]) {
                fprintf(fp1, ", %.3f", ct[k]);
            } else {
                fprintf(fp1, ", ");
            }
        }
    } else {
        fprintf(fp1, ",,,,,"); 
    }

    if ((trop_ctype&0x02) == 0x02) {
        fprintf(fp1, ", %d", sz_idx_t);
        fprintf(fp1, ", %.2f", trop_ofst);
    } else {
        fprintf(fp1, ",,"); 
    }
    if (stec_ctype != 0) {
        fprintf(fp1, ", 0x%lx", cssr->net_svmask[inet]);
    }

    for (j = s = 0; j < nsat; ++j) {
        if (stec_ctype != 0) {
            if (!((cssr->net_svmask[inet] >> (nsat - 1 - j)) & 1)) {
                continue;
            }
            gnss = sys2gnss(satsys(sat[j], &prn), NULL);
            if (s != 0) {
                fprintf(fp1, "\n,,,,,,,,,,,,,,,,,,, %d, %d", gnss, prn);
            } else {
                fprintf(fp1, ", %d, %d", gnss, prn);
            }
        }
        if ((stec_ctype&0x01) == 0x01) {
            fprintf(fp1, ", 0x%02x", stec_qual[j]);
            fprintf(fp1, ", %d", stec_type_[j]);
            
            for (k=0; k<6; k++) {
                if (k < ci_num[stec_type_[j]]) {
                    if (k > 3) {
                        fprintf(fp1, ", %.3f", ci_[j][k]);
                    } else {
                        fprintf(fp1, ", %.2f", ci_[j][k]);
                    }
                } else {
                    fprintf(fp1, ", ");
                }
            }
        }
        if (stec_ctype != 0) {
            s++;
        }
        if ((stec_ctype&0x02) == 0x02) {
            fprintf(fp1, ", %d", sz_idx_s[j]);
        }
    }
    if (stec_ctype == 0) {
        fprintf(fp1, ",,,,,,,,,,,,");
    }

    fprintf(fp2, ", 0x%lx", cssr->net_svmask[inet]);
    fprintf(fp2, ", %d", ssrg->ngp);

    for (j = 0; j < ssrg->ngp; ++j) {
        if (j != 0) {
            fprintf(fp2, "\n,,,,,,,,, %d", j + 1);
        } else {
            fprintf(fp2, ", %d", j + 1);
        }

        if (trop_ctype != 0) {
            if (total[j] != INVALID_VALUE) {
                fprintf(fp2, ", %.3f", total[j]-wet[j]);
            } else {
                fprintf(fp2, ", #N/A");
            }
            if (wet[j] != INVALID_VALUE) {
                fprintf(fp2, ", %.3f", wet[j]);
            } else {
                fprintf(fp2, ", #N/A");
            }
        } else {
            fprintf(fp2, ", ,");
        }

        if (stec_ctype != 0) {
            for (k = s = 0; k < nsat; ++k) {
                if (!((cssr->net_svmask[inet]  >> (nsat - 1 - k)) & 1)) {
                    continue;
                }

                gnss = sys2gnss(satsys(sat[k], &prn), NULL);
                if (s != 0) {
                    fprintf(fp2, "\n,,,,,,,,,,,, %d, %d", gnss, prn);
                } else {
                    fprintf(fp2, ", %d, %d", gnss, prn);
                }

                if (stec[j][s] != INVALID_VALUE) {
                    fprintf(fp2, ", %.4f", stec[j][s]);
                } else {
                    fprintf(fp2, ", #N/A");
                }
                ++s;
            }
        } else {
            fprintf(fp2, ",,,");
        }
    }

    fprintf(fp1, "\n");
    fprintf(fp2, "\n");
}

static int decode_cssr_atmos(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp1, FILE *fp2)
{
    int i, j, k, s, sync, tow, iod, ngnss, sat[CSSR_MAX_SV], nsat, sz_idx, sz_idx_t=0, sz_idx_s[CSSR_MAX_SV], sz;
    int trop_avail, stec_avail, trop_type=-1, stec_type=-1, inet, a=0, b=0, quality, trop_qual=-1, stec_qual[CSSR_MAX_SV], stec_type_[CSSR_MAX_SV];
    double total[CSSR_MAX_GP], wet[CSSR_MAX_GP], stec[CSSR_MAX_GP][CSSR_MAX_SV];
    double udint, stec0, ct[6]={0}, ci[6]={0}, ci_[CSSR_MAX_SV][6]={0}, trop_ofst=0.0, trop_residual, dlat, dlon, dstec;
    nav_t *nav = &rtcm->nav;
    ssrgp_t *ssrg;
    const double dstec_lsb_t[4] = {0.04,0.12,0.16,0.24};
    const int dstec_sz_t[4] = {4,4,5,7};

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp1);
    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp2);
    check_week_ref(rtcm, tow, ref_atmospheric);
    rtcm->time = gpst2time(rtcm->week_ref[ref_atmospheric], tow);
    nsat = svmask2sat(cssr->svmask, sat);
    
    trop_avail = getbitu(rtcm->buff, i, 2); i += 2;  /* troposphere correction availability */
    stec_avail = getbitu(rtcm->buff, i, 2); i += 2;  /* stec correction availability */
    inet = getbitu(rtcm->buff, i, 5); i += 5;       /* network id */
    
    trace(2, "decode_cssr_atmos:facility=%d tow=%d iod=%d net=%2d\n", l6facility[chidx]+1, tow, iod, inet);
    if (cssr->l6facility != l6facility[chidx]) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility[chidx]);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }
    
    ssrg = &rtcm->ssrg[inet];
    
    ssrg->ngp = getbitu(rtcm->buff, i, 6); i += 6;
    ssrg->t0 = rtcm->time;
    ssrg->network = inet;
    
    for (j = 0; j < RTCM_SSR_MAX_GP; ++j) {
        ssrg->trop_total[j] = INVALID_VALUE;
        ssrg->trop_wet[j] = INVALID_VALUE;
        ssrg->gp[j].pos[0] = 0.0;
        ssrg->gp[j].pos[1] = 0.0;
        ssrg->gp[j].pos[2] = 0.0;
        ssrg->gp[j].network = 0;
        ssrg->gp[j].update = 0;
        ssrg->nsv[j] = 0;
    }
    
    if (trop_avail != 0) {
        a = getbitu(rtcm->buff, i, 3); i += 3;
        b = getbitu(rtcm->buff, i, 3); i += 3;
        ssrg->quality = decode_cssr_quality_trop(a, b);
        trop_qual = (a<<3)|b;
    }

    if ((trop_avail&0x01) == 0x01) {
        trop_type = getbitu(rtcm->buff, i, 2); i += 2;
        for (k=0;k<4;k++) ct[k] = 0.0;
        ct[0] = decode_sval(rtcm->buff, i, 9, 0.004); i+=9;
        if (trop_type>0) {
            ct[1] = decode_sval(rtcm->buff, i, 7, 0.002); i+=7;
            ct[2] = decode_sval(rtcm->buff, i, 7, 0.002); i+=7;
        }
        if (trop_type>1) {
            ct[3] = decode_sval(rtcm->buff, i, 7, 0.001); i+=7;
        }
    }
    
    for (j = 0; j < ssrg->ngp; ++j) {
        ssrg->gp[j].pos[0] = clas_grid[inet][j][0] * D2R;
        ssrg->gp[j].pos[1] = clas_grid[inet][j][1] * D2R;
        ssrg->gp[j].pos[2] = clas_grid[inet][j][2];
        
        dlat = (ssrg->gp[j].pos[0] - ssrg->gp[0].pos[0]) * R2D;
        dlon = (ssrg->gp[j].pos[1] - ssrg->gp[0].pos[1]) * R2D;
        
        ssrg->gp[j].network = inet;
        ssrg->gp[j].update = 1;
        
        ssrg->trop_total[j] = CSSR_TROP_HS_REF + ct[0];
        if (trop_type > 0) {
            ssrg->trop_total[j] += (ct[1] * dlat) + (ct[2] * dlon);
        }
        if (trop_type > 1) {
            ssrg->trop_total[j] += ct[3] * dlat * dlon;
        }
    }

    if ((trop_avail&0x02) == 0x02) {
        sz_idx_t = getbitu(rtcm->buff, i, 1); i += 1;
        trop_ofst = getbitu(rtcm->buff, i, 4) * 0.02; i += 4;
        trace(3, "decode_cssr_atmos: network=%d, tow=%d, trop=0x%02x, stec=0x%02x, ngp=%d, trop_type=%d, ct=%.3f %.3f %.3f %.3f, sz_idx_t=%d, offset=%.3f\n",
            ssrg->network, tow, trop_avail, stec_avail, ssrg->ngp, trop_type, ct[0], ct[1], ct[2], ct[3], sz_idx_t, trop_ofst);
        sz = (sz_idx_t==0) ? 6:8;
        
        for (j = 0; j < ssrg->ngp; ++j) {
            trop_residual = decode_sval(rtcm->buff, i, sz, 0.004); i+=sz;
            if (trop_residual != INVALID_VALUE) {
                ssrg->trop_wet[j] = trop_residual + trop_ofst; 
                ssrg->trop_total[j] += ssrg->trop_wet[j];
                total[j] = ssrg->trop_total[j];
                wet[j] = ssrg->trop_wet[j];
            } else {
                ssrg->trop_total[j] = INVALID_VALUE;
                total[j] = INVALID_VALUE;
                wet[j] = INVALID_VALUE;
                trace(2,"trop(wet) is invalid: tow=%d, inet=%d, grid=%d\n",tow,inet,j);
            }
            trace(3, "decode_cssr_atmos: pos=%.3f %.3f %.3f, total=%.3f, wet=%.3f\n", ssrg->gp[j].pos[0] * R2D,
                ssrg->gp[j].pos[1] * R2D, ssrg->gp[j].pos[2], ssrg->trop_total[j], ssrg->trop_wet[j]);
        }

    }
    
    if (stec_avail != 0) {
        cssr->net_svmask[inet] = getbitu(rtcm->buff, i, nsat); i += nsat; /* stec correction type */
        trace(4, "decode_cssr_atmos: mask=0x%x\n", cssr->net_svmask[inet]);
    }
    for (j = s = 0; j < nsat; ++j) {
        if (stec_avail != 0) {
            if (!((cssr->net_svmask[inet] >> (nsat - 1 - j)) & 1)) {
                continue;
            }
            a = getbitu(rtcm->buff, i, 3); i+=3;
            b = getbitu(rtcm->buff, i, 3); i+=3;
            quality = decode_cssr_quality_stec(a,b);
        }
        for (k=0;k<6;k++) ci[k]=ci_[j][k]=0.0;
        if ((stec_avail&0x01) == 0x01) {
            stec_qual[j] = (a<<3)|b;
            stec_type_[j] = stec_type = getbitu(rtcm->buff, i, 2); i += 2;
            
            for (k=0;k<6;k++) ci[k]=0.0;
            ci_[j][0] = ci[0] = decode_sval(rtcm->buff, i, 14, 0.05); i+=14;
            if (stec_type>0) {
                ci_[j][1] = ci[1] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
                ci_[j][2] = ci[2] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
            }
            if (stec_type>1) {
                ci_[j][3] = ci[3] = decode_sval(rtcm->buff, i, 10, 0.02); i+=10;
            }
            if (stec_type>2) {
                ci_[j][4] = ci[4] = decode_sval(rtcm->buff, i, 8, 0.005); i+=8;
                ci_[j][5] = ci[5] = decode_sval(rtcm->buff, i, 8, 0.005); i+=8;
            }
        }
        if (stec_avail != 0) {
            if ((stec_avail&0x02) == 0x02) {
                sz_idx_s[j] = sz_idx = getbitu(rtcm->buff,i,2); i+=2;
                trace(3, "decode_cssr_atmos: stec_type=%d, ct=%.2f %.2f %.2f %.2f %.2f %.2f, sz_idx=%d\n",
                    stec_type, ci[0], ci[1], ci[2], ci[3], ci[4], ci[5], sz_idx);
            }
            
            for (k = 0; k < ssrg->ngp; ++k) {
                dlat = (ssrg->gp[k].pos[0] - ssrg->gp[0].pos[0]) * R2D;
                dlon = (ssrg->gp[k].pos[1] - ssrg->gp[0].pos[1]) * R2D;
                
                dstec=0.0;
                if ((stec_avail&0x02) == 0x02) {
                    dstec = decode_sval(rtcm->buff, i, dstec_sz_t[sz_idx], dstec_lsb_t[sz_idx]);
                    i += dstec_sz_t[sz_idx];
                    trace(5,"sz_idx=%d, dstec=%f\n",dstec_sz_t[sz_idx],dstec);
                }

                stec0 = ci[0];
                if (stec_type > 0) {
                    stec0 += (ci[1] * dlat) + (ci[2] * dlon);
                }
                if (stec_type > 1) {
                    stec0 += ci[3] * dlat * dlon;
                }
                if (stec_type > 2) {
                    stec0 += (ci[4] * dlat * dlat) + (ci[5] * dlon * dlon);
                }
                if (dstec == INVALID_VALUE) {
                    trace(2, "dstec is invalid: tow=%d, inet=%d, grid=%d, sat=%d, dstec=%d\n",
                        tow, inet, k, sat[j], (int)dstec);
                    ssrg->stec[k][s] = INVALID_VALUE;
                    stec[k][s] = INVALID_VALUE;
                } else {
                    ssrg->stec[k][s] = stec0 + dstec;
                    stec[k][s] = stec0 + dstec;
                }
                ssrg->stec0[k][s] = stec0;
                ssrg->sat[k][s] = sat[j];
                ++ssrg->nsv[k];
                trace(3, "decode_cssr_atmos: sat=%d, grid=%d, stec=%.4f\n",
                    ssrg->sat[k][s], k + 1, ssrg->stec[k][s]);
            }
            s++;
        }
    }

    if (stec_avail != 0) {
        for (k = 0; k < ssrg->ngp; ++k) {
            trace(4, "decode_cssr_atmos: grid=%d, nsv=%d\n", k + 1, ssrg->nsv[k]);
        }
    }

    output_cssr_atmos(rtcm, cssr, trop_avail, stec_avail, inet, trop_qual, trop_type, sz_idx_t, 
                      sz_idx_s, trop_ofst, ct, ci_, stec_qual, stec_type_, total, wet, stec, fp1, fp2);
    ssrg->update = 1;
    nav->updateac=1;
    
    check_cssr_changed_facility(cssr->l6facility);
    set_cssr_latest_trop(ssrg->t0, ssrg, inet);
    set_cssr_bank_trop(ssrg->t0, ssrg, inet);
    
    rtcm->nbit = i;
    trace(3, "decode_cssr_atmos(): tow=%d, net=%d, bits=%d\n", tow, inet, i - i0);
    return sync ? 0: 10;
}

static int check_bit_width_atmos(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int trop_avail,stec_avail,trop_type,stec_type,ngp,sz_idx,sz,j,nsat;
    uint64_t net_svmask;
    const int dstec_sz_t[4] = {4,4,5,7};
    const int trop_sz_t[3] = {9,23,30};
    const int stec_sz_t[4] = {14,38,48,64};
    
    nsat = svmask2sat(cssr->svmask,NULL);
    if (i0+36>rtcm->havebit) return FALSE;
    i0+=21;

    trop_avail = getbitu(rtcm->buff, i0, 2); i0+=2;
    stec_avail = getbitu(rtcm->buff, i0, 2); i0+=2;
    i0 += 5;
    ngp = getbitu(rtcm->buff, i0, 6); i0+=6;
    
    if (trop_avail != 0) {
        if (i0+6>rtcm->havebit) return FALSE;
        i0+=6;
    }
    if ((trop_avail&0x01) == 0x01) {
        if (i0+2>rtcm->havebit) return FALSE;
        trop_type = getbitu(rtcm->buff, i0, 2); i0+=2;
        sz = trop_sz_t[trop_type];
        if (i0+sz>rtcm->havebit) return FALSE;
        i0+=sz;
    }
    if ((trop_avail&0x02) == 0x02) {
        if (i0+5>rtcm->havebit) return FALSE;
        sz_idx = getbitu(rtcm->buff, i0, 1); i0+=1;
        i0+=4;
        sz = (sz_idx==0)?6:8;
        if (i0+sz*ngp>rtcm->havebit) return FALSE;
        i0+=sz*ngp;
    }
    
    if (stec_avail != 0) {
        if (i0+nsat>rtcm->havebit) return FALSE;
        net_svmask = getbitu(rtcm->buff, i0, nsat); i0+=nsat;
        for (j=0;j<nsat;j++) {
            if (!((net_svmask>>(nsat-1-j))&1)) continue;
            if (i0+6>rtcm->havebit) return FALSE;
            i0+=6;
            if ((stec_avail&0x01) == 0x01) {
                if (i0+2>rtcm->havebit) return FALSE;
                stec_type = getbitu(rtcm->buff, i0, 2); i0+=2;
                sz = stec_sz_t[stec_type];
                if (i0+sz>rtcm->havebit) return FALSE;
                i0+=sz;
            }
            if ((stec_avail&0x02) == 0x02) {
                if (i0+2>rtcm->havebit) return FALSE;
                sz_idx = getbitu(rtcm->buff, i0, 2); i0+=2;
                i0+=ngp*dstec_sz_t[sz_idx];
                if (i0>rtcm->havebit) return FALSE;
            }
        }
    }
    trace(4,"check_bit_width_atmos(): i0=%d, havebit=%d\n",i0,rtcm->havebit);
    return TRUE;
}

/*
 * decode service information message
 */
static int decode_cssr_si(rtcm_t *rtcm, cssr_t *cssr, int i0, int header)
{
    int i,j,sync;

    i=i0;
    sync = getbitu(rtcm->buff, i, 1); i+=1; /* multiple message indicator */
    cssr->si_cnt = getbitu(rtcm->buff, i, 3); i+=3;  /* information message counter */
    cssr->si_sz = getbitu(rtcm->buff, i, 2); i+=2; /* data size */

    for (j=0;j<cssr->si_sz;j++) {
        cssr->si_data[j] = (uint64_t)getbitu(rtcm->buff, i, 8)<<32; i+=8;
        cssr->si_data[j] |= getbitu(rtcm->buff, i, 32); i+=32;
    }
    rtcm->nbit = i;

    return sync ? 0: 10;
}

/* check if the buffer length is sufficient to decode the service information message */
static int check_bit_width_si(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int data_sz=0;
    if (i0+6>rtcm->havebit) return FALSE;
    i0+=4;
    data_sz = getbitu(rtcm->buff, i0, 2); i0+=2;

    return i0+40*(data_sz+1)<=rtcm->havebit;
}

/* decode type 4073: Melco proprietary messages */
int decode_cssr(rtcm_t *rtcm, int head)
{
    int i=12, ret = 0;
    static cssr_t _cssr = {0,};
    cssr_t *cssr = &_cssr;

    i += (head) ? 24:0;
    rtcm->subtype = getbitu(rtcm->buff,i,4); i+= 4;
    trace(4,"decode_cssr subtype=%d\n",rtcm->subtype);

    switch (rtcm->subtype) {
    case CSSR_TYPE_MASK:
        ret=decode_cssr_mask(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_OC:
        ret=decode_cssr_oc(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_CC:
        ret=decode_cssr_cc(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_CB:
        ret=decode_cssr_cb(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_PB:
        ret=decode_cssr_pb(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_BIAS:
        ret=decode_cssr_bias(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_URA:
        ret=decode_cssr_ura(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_STEC:
        ret=decode_cssr_stec(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_GRID:
        ret=decode_cssr_grid(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_COMBO:
        ret=decode_cssr_combo(rtcm, cssr, i, head, NULL);
        break;
    case CSSR_TYPE_ATMOS:
        ret = decode_cssr_atmos(rtcm, cssr, i, head, NULL, NULL);
        break;
    case CSSR_TYPE_SI:
        ret = decode_cssr_si(rtcm, cssr, i, head);
        break;
    default: break;
    }

    return ret;
}

/* read and decode QZS L6 stream */
int read_qzs_msg(rtcm_t *rtcm, unsigned char *pbuff, int nframe)
{
    int i=0, j, k, jn, jn0 = -1, prn, msgid, alert;
    unsigned char *buff;

    for (j=0;j<5;j++) {
        buff = pbuff+j*BLEN_MSG;

        prn = buff[4];
        msgid = buff[5];
        alert = (buff[6]>>7) & 0x1;

        if (msgid & 0x1) { /* head */
            jn0 = j;
            break;
        }
    }

    if (jn0 < 0) {
        memset(rtcm->buff, 0x00, sizeof(rtcm->buff));
        return 0;
    }

    trace(4,"decode_cssr prn=%d msgid=%d alert=%d\n",prn,msgid,alert);

    for (j=0,jn=jn0;j<5;j++,jn++) {
        if (jn >= 5) {
            jn = 0;
        }
        buff = pbuff+jn*BLEN_MSG;

        setbitu(rtcm->buff,i,7, buff[6] & 0x7f); i+=7;
        for (k=0;k<211;k++) {
            setbitu(rtcm->buff,i,8, buff[7+k]); i+=8;
        }
    }

    return 0;
}

/* read list of grid position from ascii file */
extern int read_grid_def(const char *gridfile)
{
    extern int gridsel;

    int no, lath, latm, lonh, lonm;
    double lat, lon, alt;
    char buff[1024], *temp, *p;
    int inet, grid[CSSR_MAX_NETWORK] = {0,}, isqzss=0, ret;
    FILE *fp;
    
    for (inet = 0; inet < CSSR_MAX_NETWORK; ++inet) {
        clas_grid[inet][0][0] = -1.0;
        clas_grid[inet][0][1] = -1.0;
        clas_grid[inet][0][2] = -1.0;
    }

    trace(2, "read_grid_def(): gridfile=%s\n", gridfile);
    fp = fopen(gridfile, "r");
    if (fp == NULL) {
        return -1;
    }
    
    while (fgets(buff, sizeof(buff), fp)) {
        if (strstr(buff, "<CSSR Grid Definition>")) {
            while (fgets(buff, sizeof(buff), fp)) {
                if ((temp = strstr(buff, "<Version>"))) {
                    p = temp + 9;
                    if ((temp = strstr(buff, "</Version>"))) {
                        *temp = '\0';
                    }
                    gridsel = atoi(p);
                    trace(2, "grid definition: version=%d\n", gridsel);
                    break;
                }
            }
            break;
        } else if (strstr(buff, "Compact Network ID    GRID No.  Latitude     Longitude   Ellipsoidal height")) {
            gridsel = 3;
            isqzss = 1;
            trace(2, "grid definition: IS attached file version%d\n", gridsel);
            break;
        } else {
            trace(1, "grid definition: invalid format%d\n", gridsel);
            fclose(fp);
            return -1;
        }
    }
    fclose(fp);

    fp = fopen(gridfile, "r");
    if (fp == NULL) {
        return -1;
    }

    if (isqzss == 0) { 
        while (fgets(buff, sizeof(buff), fp)) {
            if (sscanf(buff, "<Network%d>", &inet)) {
                while (fscanf(fp, "%d\t%d\t%d\t%lf\t%d\t%d\t%lf\t%lf",
                              &no, &lath, &latm, &lat, &lonh, &lonm, &lon, &alt) > 0) {
                    if (inet >= 0 && inet < CSSR_MAX_NETWORK) {
                        clas_grid[inet][grid[inet]][0] = (double)lath + ((double)latm/60.0) + (lat/3600.0);
                        clas_grid[inet][grid[inet]][1] = (double)lonh + ((double)lonm/60.0) + (lon/3600.0);
                        clas_grid[inet][grid[inet]][2] = alt;
                        ++grid[inet];
                        clas_grid[inet][grid[inet]][0] = -1.0;
                        clas_grid[inet][grid[inet]][1] = -1.0;
                        clas_grid[inet][grid[inet]][2] = -1.0;
                    }
                }
            }
        }
    } else {
        if (!fgets(buff, sizeof(buff), fp)) return -1;
        while ( (ret=fscanf(fp, "%d %d %lf %lf %lf", &inet, &no, &lat, &lon, &alt)) != EOF ) {
            if (inet >= 0 && inet < CSSR_MAX_NETWORK && ret == 5) {
                clas_grid[inet][grid[inet]][0] = lat;
                clas_grid[inet][grid[inet]][1] = lon;
                clas_grid[inet][grid[inet]][2] = alt;
                ++grid[inet];
                clas_grid[inet][grid[inet]][0] = -1.0;
                clas_grid[inet][grid[inet]][1] = -1.0;
                clas_grid[inet][grid[inet]][2] = -1.0;
            }
            trace(3, "grid_info(fscanf:%d), %d, %d, %lf, %lf, %lf\n", ret, inet, no, lat, lon, alt);
        }
    }
    fclose(fp);
    return 0;
}

/* decode QZS L6 CLAS stream */
extern int decode_qzs_msg(rtcm_t *rtcm, int head, uint8_t *frame, FILE **ofp)
{
    static int startdecode[L6_CH_NUM] = {0,};
    static cssr_t _cssr[L6_CH_NUM] = {0,};

    cssr_t *cssr = &_cssr[chidx];
    int startbit, week;
    int i, ret = 0;
    double tow;

    if (*frame == 0x00 || rtcm->nbit == -1) {
        return 0;
    }
    
    i = startbit = rtcm->nbit;
    if ((i + 16) > rtcm->havebit) {
        return 0;
    }
    rtcm->ctype = getbitu(rtcm->buff, i, 12); i+= 12;
    if (rtcm->ctype != 4073) {
        trace(4, "cssr: decode terminate: frame=%02x, havebit=%d, nbit=%d\n", *frame, rtcm->havebit, rtcm->nbit);
        rtcm->nbit = -1;
        *frame = 0;
        return 0;
    }
    rtcm->subtype = getbitu(rtcm->buff, i, 4); i+= 4;
    if (rtcm->subtype != 1 && startdecode[chidx] == FALSE) {
        return 0;
    }

    switch (rtcm->subtype) {
    case CSSR_TYPE_MASK:
        if (!check_bit_width_mask(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_OC:
        if (!check_bit_width_oc(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_CC:
        if (!check_bit_width_cc(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_CB:
        if (!check_bit_width_cb(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_PB:
        if (!check_bit_width_pb(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_BIAS:
        if (!check_bit_width_bias(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_URA:
        if (!check_bit_width_ura(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_STEC:
        if (!check_bit_width_stec(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_GRID:
        if (!check_bit_width_grid(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_COMBO:
        if (!check_bit_width_combo(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_ATMOS:
        if (!check_bit_width_atmos(rtcm, cssr, i)) return FALSE;
        break;
    case CSSR_TYPE_SI:
        if (!check_bit_width_si(rtcm, cssr, i)) return FALSE;
        break;
    case 0:
        trace(1, "invalid process: frame=%02x, havebit=%d, nbit=%d\n", *frame, rtcm->havebit, rtcm->nbit);
        return 0;
    }
    
    trace(4, "cssr: frame=0x%02x, ctype=%d, subtype=%d\n", *frame, rtcm->ctype, rtcm->subtype);
    tow = time2gpst(timeget(), &week);
    
    switch (rtcm->subtype) {
        case CSSR_TYPE_MASK:
            ret=decode_cssr_mask(rtcm, cssr, i, head, ofp[0]);
            if (startdecode[chidx] == FALSE) {
                trace(1, "start CSSR decoding: week=%d, tow=%.1f\n", week, tow);
                startdecode[chidx] = TRUE;
            }
            break;
        case CSSR_TYPE_OC:
            ret=decode_cssr_oc(rtcm, cssr, i, head, ofp[1]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_CC:
            ret=decode_cssr_cc(rtcm, cssr, i, head, ofp[2]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_CB:
            ret=decode_cssr_cb(rtcm, cssr, i, head, ofp[3]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_PB:
            ret=decode_cssr_pb(rtcm, cssr, i, head, ofp[4]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_BIAS:
            ret=decode_cssr_bias(rtcm, cssr, i, head, ofp[5]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_URA:
            ret=decode_cssr_ura(rtcm, cssr, i, head, ofp[6]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_STEC:
            ret=decode_cssr_stec(rtcm, cssr, i, head, ofp[7]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_GRID:
            ret=decode_cssr_grid(rtcm, cssr, i, head, ofp[8]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_COMBO:
            ret=decode_cssr_combo(rtcm, cssr, i, head, ofp[11]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_ATMOS:
            ret=decode_cssr_atmos(rtcm, cssr, i, head, ofp[12], ofp[13]);
            if (ret == -1) {
                rtcm->nbit = -1;
                cssr->iod = -1;
                *frame = 0;
                return 0;
            }
            break;
        case CSSR_TYPE_SI:
            ret=decode_cssr_si(rtcm, cssr, i, head);
            break;
        default: break;
    }

    return ret;
}

static int GetCSSRTime(unsigned char *buff, int *value)
{
    static int basetow = -1;
    static int savetow = -1;
    int i = 32 + 8 + 8 + 1;
    int type, subtype, tow;

    type = getbitu(buff, i, 12); i += 12;
    *value = -1;
    if (type != 4073) {
        return -1;
    }
    subtype = getbitu(buff, i, 4); i += 4;

    if (subtype ==  CSSR_TYPE_MASK) {
        tow = getbitu(buff, i, 20); i += 20;
        basetow = (int)floor((double)tow / 3600.0) * 3600;
        *value = subtype;
    } else {
        tow = getbitu(buff, i, 12); i +=12;
        if (basetow != -1) {
            tow += basetow;
            if ((tow < savetow) && (savetow - tow) < 3600*24*7/2) {
                if ((tow < savetow) && (savetow - tow) > 3000) {
                    return -1;
                }
            }
            savetow = tow;
        }
    }
    return tow;
}

static void output_cssr_header(int chidx, int *nframe, unsigned char buff[][BLEN_MSG], int *tow, FILE *fp) 
{
    char str1[256] = {'\0'}, cVenderid[256] = {'\0'}, cMsggenid[256] = {'\0'}, cTransptn[256] = {'\0'}, cSubframe[256] = {'\0'};
    int msgid, alert, binary, base, venderid, msggenid, transptn, subframe;
    int type;
    char *endp;

    if (fp == NULL) return;

    if (nframe[chidx] == 0) {
        tow[chidx] = GetCSSRTime(buff[chidx], &type);
    }
    
    if (tow[chidx] != -1) {
        fprintf(fp, "%d, 0x%02x-0x%02x-0x%02x-0x%02x, ", tow[chidx] + nframe[chidx], buff[chidx][0], buff[chidx][1], buff[chidx][2], buff[chidx][3]);
    } else {
        fprintf(fp, "#N/A, 0x%02x-0x%02x-0x%02x-0x%02x, ", buff[chidx][0], buff[chidx][1], buff[chidx][2], buff[chidx][3]);
    }
    binary = 0;
    base = 1; 
    msgid = buff[chidx][5];
    while(msgid>0){ 
        binary = binary + ( msgid % 2 ) * base; 
        msgid = msgid / 2; 
        base = base * 10;
    }
    sprintf(str1, "%d", binary);
    memcpy(cVenderid, &str1[0], 3);
    memcpy(cMsggenid, &str1[3], 2);
    memcpy(cTransptn, &str1[5], 2);
    memcpy(cSubframe, &str1[7], 1);

    venderid = strtol(cVenderid, &endp, 2);
    msggenid = strtol(cMsggenid, &endp, 2);
    transptn = strtol(cTransptn, &endp, 2);
    subframe = strtol(cSubframe, &endp, 2);
    fprintf(fp, "%d, %d, %d, %d, %d, %d, %d\n", buff[chidx][4], buff[chidx][5], venderid, msggenid, transptn, subframe, alert = (buff[chidx][6]>>7) & 0x1);
    
    return;
}

extern int get_correct_fac(int msgid) {
    int fac = (msgid & 0x18) >> 3;
    int ptn = (msgid & 0x06) >> 1;

    if (ptn == 0) {
        switch(fac) {
            case 0: return 0;
            case 1: return 2;
            case 2: return 1;
            case 3: return 3;
            default: return -1;
        }
    } else {
        switch(fac) {
            case 0: return 2;
            case 1: return 0;
            case 2: return 3;
            case 3: return 1;
            default: return -1;
        }
    }
}

/* decode cssr messages in the QZS L6 subframe */
extern int input_cssr(rtcm_t *cssr, unsigned char data, uint8_t *frame, FILE **ofp)
{
    static uint32_t preamble[L6_CH_NUM] = {0,};
    static uint64_t data_p[L6_CH_NUM] = {0,};
    uint8_t prn, msgid, alert;
    static int nframe[L6_CH_NUM] = {0,};
    static unsigned char buff[L6_CH_NUM][BLEN_MSG];
    static int decode_start[L6_CH_NUM] = {0,};
    static int tow[L6_CH_NUM];

    trace(5,"input_cssr: data=%02x\n",data);

    /* synchronize frame */
    if (cssr->nbyte==0) {
        preamble[chidx] = (preamble[chidx] << 8) | data;
        data_p[chidx] = (data_p[chidx] << 8) | data;
        if (preamble[chidx] != L6FRMPREAMB) {
            return 0;
        }
        preamble[chidx] = 0;
        buff[chidx][cssr->nbyte++]=(L6FRMPREAMB>>24) & 0xff;
        buff[chidx][cssr->nbyte++]=(L6FRMPREAMB>>16) & 0xff;
        buff[chidx][cssr->nbyte++]=(L6FRMPREAMB>>8) & 0xff;
        buff[chidx][cssr->nbyte++]=data;
        return 0;
    }
    buff[chidx][cssr->nbyte++]=data;
    cssr->len = BLEN_MSG;

    if (cssr->nbyte<cssr->len) return 0;
    cssr->nbyte=0;

    prn = buff[chidx][4];
    msgid = buff[chidx][5];
    alert = (buff[chidx][6]>>7) & 0x1;
    if (alert != 0) {
        trace(1, "CSSR frame alert!: tow=%.1f\n", time2gpst(timeget(), NULL));
    }

    l6delivery[chidx] = prn;
    l6facility[chidx] = get_correct_fac(msgid);

    if (msgid & 0x01) { /* Subframe indicator */
        if (decode_start[chidx] == 0) {
            trace(1, "CSSR frame first recieve: tow=%.1f\n", time2gpst(timeget(), NULL));
        }
        cssr->havebit = 0;
        decode_start[chidx] = 1;
        cssr->nbit = 0;
        *frame = 0;
        nframe[chidx] = 0;
    } else if (nframe[chidx] >= 5) {
        return 0;
    }
    if (decode_start[chidx] == 1) {
        int i=1695*nframe[chidx],j;

        output_cssr_header(chidx,nframe,buff,tow,ofp[10]);        
        
        setbitu(cssr->buff, i, 7, buff[chidx][6] & 0x7f); i+=7;
        for (j = 0; j < 211; j++) {
            setbitu(cssr->buff, i, 8, buff[chidx][7+j]); i+=8;
        }
        cssr->havebit += 1695;
        *frame |= (1<<nframe[chidx]);
        nframe[chidx]++;
    }
    return 0;
}


/* decode cssr messages from file stream ---------------------------------------------*/
extern int input_cssrf(rtcm_t *cssr, FILE *fp, FILE **ofp)
{
    static uint8_t frame[L6_CH_NUM] = {0,};
    int i,data=0,ret;

    trace(4,"input_cssrf: data=%02x\n",data);

    for (i=0;i<4096;i++) {
        if ((ret=decode_qzs_msg(cssr, 0, &frame[chidx], ofp))) return ret;
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_cssr(cssr, (unsigned char)data, &frame[chidx],ofp))) return ret;
    }
    return 0; /* return at every 4k bytes */
}

/* open QZSS L6 message file -------------------------------------------------*/
extern FILE *open_L6(char **infile, int n)
{
    FILE *fp;
    char *ext;
    int i;
    
    for (i = 0; i < n; i++) {
        if (!(ext = strrchr(infile[i], '.'))) continue;
        if (!strcmp(ext, ".l6") || !strcmp(ext, ".L6")) {
            break;
        }
    }
    if (i >= n) {
        fprintf(stderr, "No L6 message file in input files.\n");
        return NULL;
    }
    if (!(fp = fopen(infile[i], "rb"))) {
        fprintf(stderr, "L6 message file open error. %s\n", infile[i]);
        return NULL;
    }
    return fp;
}

extern int dumpcssr(char **infile, int n, FILE **ofp, filopt_t *fopt){

    int ret;
    static rtcm_t rtcm;
    FILE *fp;

    init_rtcm(&rtcm);

    /* open QZSS L6 message file */
    if (!(fp = open_L6(infile, n))) {
        free_rtcm(&rtcm);
        return 0;
    }
    /* read grid definition file */
    if (read_grid_def(fopt->grid)) {
        fprintf(stderr, "Grid file read error. %s\n", fopt->grid);
        showmsg("Grid file read error. %s\n", fopt->grid);
        fclose(fp);
        return -1;
    }

    while (1) {
        if ((ret = input_cssrf(&rtcm, fp, ofp)) < -1) {
            break;
        }
    }
    fclose(fp);
    return 0;
}

extern int open_outputfiles(FILE **ofp)
{
    char filename[512];
    static char parsename[512] = "parse_cssr";

    /* mask */
    sprintf(filename, "%s_type1.csv", parsename);
    if (!(ofp[0] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[0], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,No. of GNSS,");
    fprintf(ofp[0], "GNSS ID,Compact SSR Satellite Mask,Compact SSR Signal Mask,Cell-Mask Availability Flag,");
    fprintf(ofp[0], "[Satellite Number],Compact SSR Cell mask\n");

    /* orbit */
    sprintf(filename, "%s_type2.csv", parsename);
    if (!(ofp[1] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[1], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[1], "[GNSS ID],[Satellite Number],GNSS IODE,Compact SSR Delta Radial,Compact SSR Delta Along-Track,");
    fprintf(ofp[1], "Compact SSR Delta Cross-Track\n");

    /* clock */
    sprintf(filename, "%s_type3.csv", parsename);
    if (!(ofp[2] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[2], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[2], "[GNSS ID],[Satellite Number],Compact SSR Delta Clock C0\n");

    /* code bias */
    sprintf(filename, "%s_type4.csv", parsename);
    if (!(ofp[3] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[3], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[3], "[GNSS ID],[Satellite Number],[Satellite Signal],Compact SSR Code Bias\n");

    /* phase bias */
    sprintf(filename, "%s_type5.csv", parsename);
    if (!(ofp[4] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[4], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[4], "[GNSS ID],[Satellite Number],[Satellite Signal],Compact SSR Phase Bias,");
    fprintf(ofp[4], "Compact SSR Phase Discontinuity Indicator\n");

    /* bias */
    sprintf(filename, "%s_type6.csv", parsename);
    if (!(ofp[5] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[5], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[5], "Code Bias Existing Flag,Phase Bias Existing Flag,Network Bias Correction,Compact Network ID,Network SV Mask,");
    fprintf(ofp[5], "[GNSS ID],[Satellite Number],[Satellite Signal],Compact SSR Code Bias,Compact SSR Phase Bias,");
    fprintf(ofp[5], "Compact SSR Phase Discontinuity Indicator\n");

    /* URA */
    sprintf(filename, "%s_type7.csv", parsename);
    if (!(ofp[6] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[6], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[6], "[GNSS ID],[Satellite Number],SSR URA\n");

    /* STEC */
    sprintf(filename, "%s_type8.csv", parsename);
    if (!(ofp[7] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[7], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[7], "Compact SSR STEC Correction Type,Compact Network ID,Network SV Mask,");
    fprintf(ofp[7], "[GNSS ID],[Satellite Number],SSR STEC Quality Indicator,Polynomial Coefficients C00,Polynomial Coefficients C01,");
    fprintf(ofp[7], "Polynomial Coefficients C10,Polynomial Coefficients C11\n");

    /* grid */
    sprintf(filename,  "%s_type9.csv", parsename);
    if (!(ofp[8] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[8], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[8], "Tropospheric Delay Correction Type,STEC Residual Correction Range,Compact Network ID,");
    fprintf(ofp[8], "Network SV Mask,Tropospheric Delay Quality Indicator,No. of Grids,[Grid Number],");
    fprintf(ofp[8], "Troposphere Hydro-Static Vertical Delay,Troposphere Wet Vertical Delay,[GNSS ID],");
    fprintf(ofp[8], "[Satellite Number],STEC Residual Correction\n");

    /* combo */
    sprintf(filename,  "%s_type11.csv", parsename);
    if (!(ofp[11] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[11], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[11], "Orbit Existing Flag,Clock Existing Flag,Network Correction,Network ID,Network SV Mask,");
    fprintf(ofp[11], "[GNSS ID],[Satellite Number],GNSS IODE,Compact SSR Delta Radial,Compact SSR Delta Along-Track,");
    fprintf(ofp[11], "Compact SSR Delta Cross-Track,Compact SSR Delta Clock C0\n");

    /* atmos */
    sprintf(filename,  "%s_type12_stec.csv", parsename);
    if (!(ofp[12] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[12], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[12], "Tropospheric Correction Availability,STEC Correction Availability,Compact Network ID,No. of Grids,");
    fprintf(ofp[12], "Troposphere Quality Indicator,Tropospheric Correction Type,Troposphere Polynomial Coefficients T00,");
    fprintf(ofp[12], "Troposphere Polynomial Coefficients T01,Troposphere Polynomial Coefficients T10,");
    fprintf(ofp[12], "Troposphere Polynomial Coefficients T11,Troposphere Residual Size,Troposphere Residual Offset,");
    fprintf(ofp[12], "Network SV Mask,[GNSS ID],[Satellite Number],STEC Quality Indicator,STEC Correction Type,");
    fprintf(ofp[12], "STEC Polynomial Coefficients C00,STEC Polynomial Coefficients C01,STEC Polynomial Coefficients C10,");
    fprintf(ofp[12], "STEC Polynomial Coefficients C11,STEC Polynomial Coefficients C02,STEC Polynomial Coefficients C20,");
    fprintf(ofp[12], "STEC Residual Size\n");

    sprintf(filename,  "%s_type12_grid.csv", parsename);
    if (!(ofp[13] = fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[13], "Message Number,Message Sub Type ID,GNSS Epoch Time 1s,SSR Update Interval,Multiple Message Indicator,IOD SSR,");
    fprintf(ofp[13], "Compact Network ID,Network SV Mask,No. of Grids,[Grid Number],Troposphere Hydro-Static Vertical Delay,");
    fprintf(ofp[13], "Troposphere Wet Vertical Delay,[GNSS ID],[Satellite Number],STEC Residual Correction[TECU]\n");

    /* header */
    sprintf(filename, "%s_header.csv", parsename);
    if (!(ofp[10]=fopen(filename, "w"))) {
        return -1;
    }
    fprintf(ofp[10], "Epoch Time,Preamble,PRN,L6 message type ID,Vender ID,Message Generation Facility ID and CLAS Transmit Pattern ID,CLAS Transmit Pattern ID,Subframe indicator,Alert Flag\n");

    return 0;
}

extern void close_outputfiles(FILE **ofp)
{
    int i;
    for (i=0; i<CSSR_TYPE_NUM; i++) {
        if (ofp[i] != NULL) {
            fclose(ofp[i]);
        }
    }
}

/* init global cssr object */
extern void init_cssr_object(int idx)
{
    chidx = idx;
    memset(&cssrObject[chidx].LatestTrop, 0x00, sizeof(cssr_latest_trop));
    memset(cssrObject[chidx].OrbitBank, 0x00, sizeof(cssrObject[chidx].OrbitBank));
    memset(cssrObject[chidx].ClockBank, 0x00, sizeof(cssrObject[chidx].ClockBank));
    memset(cssrObject[chidx].BiasBank, 0x00, sizeof(cssrObject[chidx].BiasBank));
    memset(cssrObject[chidx].TropBank, 0x00, sizeof(cssrObject[chidx].TropBank));
    cssrObject[chidx].Facility = -1;
    cssrObject[chidx].separation = 0;
    cssrObject[chidx].NextOrbit = 0;
    cssrObject[chidx].NextClock = 0;
    cssrObject[chidx].NextBias = 0;
    cssrObject[chidx].NextTrop = 0;
    cssrObject[chidx].use = 1;
}

/* set current cssr object index */
extern void set_cssr_ch_idx(int ch)
{
    chidx = ch;
}
/* get current cssr object index */
extern int get_cssr_ch_idx(void)
{
    return chidx;
}
