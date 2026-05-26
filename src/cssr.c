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

static int l6delivery = -1;
static int l6facility = -1;

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
    if (rtcm->tow0 != -1) {
        if (rtcm->tow_ref[i] != -1 && ((tow - rtcm->tow_ref[i]) < (-86400*7/2))) {
            ++rtcm->week_ref[i];
        }
        rtcm->tow_ref[i] = tow;
    }
}

typedef struct _CSSRBank {
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
    int sisadjust;
    int facility;
    int gridnum;
    int network;
    int use;
} CSSRBank;

typedef struct _CSSROrbitBank {
    int use;
    gtime_t time;
    int validnetwork;
    int prn[CSSR_MAX_NETWORK][MAXSAT];
    double udi[MAXSAT];
    int iod[MAXSAT];
    int iode[CSSR_MAX_NETWORK][MAXSAT];
    double deph0[CSSR_MAX_NETWORK][MAXSAT];
    double deph1[CSSR_MAX_NETWORK][MAXSAT];
    double deph2[CSSR_MAX_NETWORK][MAXSAT];
} CSSROrbitBank;

typedef struct _CSSRBiasBank {
    int use;
    gtime_t time;
    int prn[CSSR_MAX_NETWORK][MAXSAT];
    double udi[MAXSAT];
    int iod[MAXSAT];
    int bflag[CSSR_MAX_NETWORK];
    int smode[CSSR_MAX_NETWORK][MAXSAT][MAXCODE];
    int sflag[CSSR_MAX_NETWORK][MAXSAT][MAXCODE];
    double cbias[CSSR_MAX_NETWORK][MAXSAT][MAXCODE];
    double pbias[CSSR_MAX_NETWORK][MAXSAT][MAXCODE];
} CSSRBiasBank;

typedef struct _CSSRClockBank {
    int use;
    gtime_t time;
    int validnetwork;
    int prn[CSSR_MAX_NETWORK][MAXSAT];
    double udi[MAXSAT];
    int iod[MAXSAT];
    gtime_t ssrtime[CSSR_MAX_NETWORK][MAXSAT];
    double c0[CSSR_MAX_NETWORK][MAXSAT];
} CSSRClockBank;

typedef struct _CSSRTropBank {
    int use;
    gtime_t time;
    int gridnum[CSSR_MAX_NETWORK];
    double gridpos[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][3];
    double total[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    double wet[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    int satnum[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    int prn[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    double iono[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
} CSSRTropBank;

typedef struct _CSSRLatestTrop {
    int gridnum[CSSR_MAX_NETWORK];
    double gridpos[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][3];
    gtime_t troptime[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    double total[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    double wet[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
    gtime_t stectime[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    int prn[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    double stec0[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
    double stec[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][MAXSAT];
} CSSRLatestTrop;

#define BIASBANKNUM  128
#define TROPBANKNUM  128
#define ORBBANKNUM   128
#define CLKBANKNUM   128

static CSSROrbitBank OrbitBank[ORBBANKNUM];
static CSSRClockBank ClockBank[CLKBANKNUM];
static CSSRBiasBank  BiasBank[BIASBANKNUM];
static CSSRTropBank  TropBank[TROPBANKNUM];
static CSSRLatestTrop LatestTrop;
static CSSRBank CurrentCSSR;
static CSSRBank BackupCSSR;
static grid_t BackupGrid;
static int Facility  = -1;
static int NextOrbit = 0;
static int NextClock = 0;
static int NextBias  = 0;
static int NextTrop  = 0;

static void CheckCSSRChangedFacility(int facility)
{
    if (Facility != facility) {
        trace(4, "CCSR bank clear: facility changed, %d ---> %d\n", Facility + 1, facility + 1);
        memset(&LatestTrop, 0x00, sizeof(CSSRLatestTrop));
        memset(OrbitBank, 0x00, sizeof(OrbitBank));
        memset(ClockBank, 0x00, sizeof(ClockBank));
        memset(BiasBank,  0x00, sizeof(BiasBank));
        memset(TropBank,  0x00, sizeof(TropBank));
        Facility = facility;
        NextOrbit = 0;
        NextClock = 0;
        NextBias  = 0;
        NextTrop  = 0;
    }
}

static CSSROrbitBank *GetSameOrbitCorrection(gtime_t time)
{
    int i;

    for (i = 0; i < ORBBANKNUM; ++i) {
        if (OrbitBank[i].use == TRUE && timediff(OrbitBank[i].time, time) == 0.0) {
            return &OrbitBank[i];
        }
    }
    return NULL;
}

static CSSROrbitBank *GetCloseOrbitCorrection(gtime_t time, int network, double age)
{
    int search = 0;
    int pos = -1;
    int i;

    for (i = 0; i < ORBBANKNUM; ++i) {
        if (OrbitBank[i].use == TRUE && (OrbitBank[i].validnetwork & (1 << network))) {
            search = network;
            break;
        }
    }
    for (i = 0; i < ORBBANKNUM; ++i) {
        if (OrbitBank[i].use == TRUE) {
            trace(5, "GetCloseOrbitCorrection(): orbit=%.1f, network=%d, diff=%.1f\n", time2gpst(OrbitBank[i].time, NULL), OrbitBank[i].validnetwork, timediff(OrbitBank[i].time, time));
        }
        if (OrbitBank[i].use == TRUE && (OrbitBank[i].validnetwork & (1 << search)) && fabs(timediff(OrbitBank[i].time, time)) <= age) {
            if (pos != -1 && fabs(timediff(OrbitBank[pos].time, time)) < fabs(timediff(OrbitBank[i].time, time))) {
                continue;
            }
            if (timediff(OrbitBank[i].time, time) > 0.0) {
                continue;
            }
            pos = i;
        }
    }
    if (pos != -1) {
        trace(4, "GetCloseOrbitCorrection(): orbit=%.1f, network=%d, diff=%.1f\n", time2gpst(OrbitBank[pos].time, NULL),
            OrbitBank[pos].validnetwork, timediff(OrbitBank[pos].time, time));
        return &OrbitBank[pos];
    }
    return NULL;
}

void SetCSSRBankOrbit(gtime_t time, nav_t *nav, int network)
{
    CSSROrbitBank *orbit;
    int i;

    if ((orbit = GetSameOrbitCorrection(time)) == NULL) {
        OrbitBank[NextOrbit].time = time;
        OrbitBank[NextOrbit].use = TRUE;
        orbit = &OrbitBank[NextOrbit];
        if (++NextOrbit >= ORBBANKNUM) {
            NextOrbit = 0;
        }
        memset(orbit->prn, 0x00, sizeof(orbit->prn));
        orbit->validnetwork = 0x00;
    }
    orbit->validnetwork |= (1 << network);

    for (i = 0; i < MAXSAT; ++i) {
        if (network == 0) {
            if (nav->ssr[i].update_oc == 1) {
                if (nav->ssr[i].deph[0] != INVALID_VALUE) {
                    trace(3, "SetCSSRBankOrbit(): tow=%.1f, prn=%2d, iode=%d, deph=%.4f\n", time2gpst(time, NULL),
                            i + 1, nav->ssr[i].iode, nav->ssr[i].deph[0]);
                    orbit->prn[network][i] = i + 1;
                    orbit->udi[i] = nav->ssr[i].udi[0];
                    orbit->iod[i] = nav->ssr[i].iod[0];
                    orbit->iode[network][i] = nav->ssr[i].iode;
                    orbit->deph0[network][i] = nav->ssr[i].deph[0];
                    orbit->deph1[network][i] = nav->ssr[i].deph[1];
                    orbit->deph2[network][i] = nav->ssr[i].deph[2];
                } else {
                    trace(3, "SetCSSRBankOrbit(): tow=%.1f, prn=%2d, iode=%d, deph=#N/A\n", time2gpst(time, NULL),
                            i + 1, nav->ssr[i].iode);
                    orbit->prn[network][i] = 0;
                }
            }
        } else {
            if (nav->extcorr[network][i].update_oc == 1) {
                if (nav->extcorr[network][i].deph[0] != INVALID_VALUE) {
                    trace(3, "SetCSSRBankOrbit(): tow=%.1f, network=%d, prn=%2d, iode=%d, deph=%.4f\n", time2gpst(time, NULL),
                            network, i + 1, nav->extcorr[network][i].iode, nav->extcorr[network][i].deph[0]);
                    orbit->prn[network][i] = i + 1;
                    orbit->udi[i] = nav->extcorr[network][i].udi[2];
                    orbit->iod[i] = nav->extcorr[network][i].iod[2];
                    orbit->iode[network][i] = nav->extcorr[network][i].iode;
                    orbit->deph0[network][i] = nav->extcorr[network][i].deph[0];
                    orbit->deph1[network][i] = nav->extcorr[network][i].deph[1];
                    orbit->deph2[network][i] = nav->extcorr[network][i].deph[2];
                } else {
                    trace(3, "SetCSSRBankOrbit(): tow=%.1f, network=%d, prn=%2d, iode=%d, deph=#N/A\n", time2gpst(time, NULL),
                            network, i + 1, nav->extcorr[network][i].iode);
                    orbit->prn[network][i] = 0;
                }
            }
        }
    }
    trace(4, "SetCSSRBankOrbit(): next=%d\n", NextOrbit);
}

static CSSRBiasBank *GetSameBiasCorrection(gtime_t time)
{
    int i;

    for (i = 0; i < BIASBANKNUM; ++i) {
        if (BiasBank[i].use == TRUE && timediff(BiasBank[i].time, time) == 0.0) {
            return &BiasBank[i];
        }
    }
    return NULL;
}

static CSSRBiasBank *GetCloseBiasCorrection2(gtime_t time, int network, double age)
{
    int pos = -1;
    int i;

    trace(3, "GetCloseBiasCorrection(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);

    if (network >= 1 && network <= CSSR_MAX_NETWORK) {
        for (i = 0; i < BIASBANKNUM; ++i) {
            if (BiasBank[i].use == TRUE) {
                trace(5, "GetCloseBiasCorrection(): bias=%.1f, diff=%.1f, network=%d, flag=0x%02x\n", time2gpst(BiasBank[i].time, NULL),
                    timediff(BiasBank[i].time, time), network, BiasBank[i].bflag[network-1]);
            }
            if (BiasBank[i].use == TRUE && timediff(time, BiasBank[i].time) >= 0.0 && timediff(time, BiasBank[i].time) <= age) {
                if (pos != -1 && timediff(BiasBank[pos].time, BiasBank[i].time) > 0.0) {
                    continue;
                }
                if (BiasBank[i].bflag[network-1] == 0x03) {
                    pos = i;
                }
            }
        }
        if (pos != -1) {
            trace(4, "GetCloseBiasCorrection(): bias=%.1f, diff=%.1f\n", time2gpst(BiasBank[pos].time, NULL),
                timediff(BiasBank[pos].time, time));
            return &BiasBank[pos];
        }
    }
    return NULL;
}

static CSSRBiasBank *GetCloseCBiasCorrection(gtime_t time, int network, double age)
{
    int pos = -1;
    int i;

    trace(3, "GetCloseCBiasCorrection(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);

    if (network >= 1 && network <= CSSR_MAX_NETWORK) {
        for (i = 0; i < BIASBANKNUM; ++i) {
            if (BiasBank[i].use == TRUE) {
                trace(5, "GetCloseCBiasCorrection(): bias=%.1f, diff=%.1f, network=%d, flag=0x%02x\n", time2gpst(BiasBank[i].time, NULL),
                    timediff(BiasBank[i].time, time), network, BiasBank[i].bflag[network-1]);
            }
            if (BiasBank[i].use == TRUE && timediff(time, BiasBank[i].time) >= 0.0 && timediff(time, BiasBank[i].time) <= age) {
                if (pos != -1 && timediff(BiasBank[pos].time, BiasBank[i].time) > 0.0) {
                    continue;
                }
                if ((BiasBank[i].bflag[network-1] & 0x01) == 0x01) {
                    pos = i;
                }
            }
        }
        if (pos != -1) {
            trace(4, "GetCloseCBiasCorrection(): bias=%.1f, diff=%.1f, flag=0x%02x\n", time2gpst(BiasBank[pos].time, NULL),
                timediff(BiasBank[pos].time, time), BiasBank[pos].bflag[network-1]);
            return &BiasBank[pos];
        }
    }
    return NULL;
}

static CSSRBiasBank *GetClosePBiasCorrection(gtime_t time, int network, double age)
{
    int pos = -1;
    int i;

    trace(3, "GetClosePBiasCorrection(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);

    if (network >= 1 && network <= CSSR_MAX_NETWORK) {
        for (i = 0; i < BIASBANKNUM; ++i) {
            if (BiasBank[i].use == TRUE) {
                trace(5, "GetClosePBiasCorrection(): bias=%.1f, diff=%.1f, network=%d, flag=0x%02x\n", time2gpst(BiasBank[i].time, NULL),
                    timediff(BiasBank[i].time, time), network, BiasBank[i].bflag[network-1]);
            }
            if (BiasBank[i].use == TRUE && timediff(time, BiasBank[i].time) <= age) {
                if (pos != -1 && timediff(BiasBank[pos].time, BiasBank[i].time) > 0.0) {
                    continue;
                }
                if ((BiasBank[i].bflag[network-1] & 0x02) == 0x02) {
                    pos = i;
                }
            }
        }
        if (pos != -1) {
            trace(4, "GetClosePBiasCorrection(): bias=%.1f, diff=%.1f, flag=0x%02x\n", time2gpst(BiasBank[pos].time, NULL),
                timediff(BiasBank[pos].time, time), BiasBank[pos].bflag[network-1]);
            return &BiasBank[pos];
        }
    }
    return NULL;
}


static int GetAddBiasPoint(CSSRBiasBank *bias, int network, int prn, int mode)
{
    int i;

    for (i = 0; i < MAXCODE; ++i) {
        if (bias->smode[network-1][prn-1][i] == mode) {
            return i;
        }
    }
    for (i = 0; i < MAXCODE; ++i) {
        if (bias->smode[network-1][prn-1][i] == 0) {
            return i;
        }
    }
    return -1;
}

double GetCSSRCBiasValue(CSSRBiasBank *bias, int network, int sat, int mode)
{
    int i;
    
    for (i = 0; i < MAXCODE; ++i) {
        if (bias->smode[network-1][sat-1][i] == mode) {
            trace(4, "GetCSSRCBiasValue(): network=%d, sat=%d, mode=%d, cbias=%f\n", network,
                bias->prn[network-1][sat-1], bias->smode[network-1][sat-1][i],
                bias->cbias[network-1][sat-1][i]);
            return bias->cbias[network-1][sat-1][i];
        }
    }
    return INVALID_VALUE;
}

int GetCSSRIODCBias(CSSRBiasBank *bias, int network)
{
    int i;
    
    for (i = 0; i < MAXSAT; ++i) {
        if (bias->prn[network-1][i] != 0) {
            trace(4, "GetCSSRIODCBias(): tow=%.1f, sat=%d, iod=%d\n",
                time2gpst(bias->time, NULL), i + 1, bias->iod[i]);
            return bias->iod[i];
        }
    }
    return -1;
}

void SetCSSRBankCBias(gtime_t time, nav_t *nav, int network, int iod)
{
    CSSRBiasBank *bias, *basebias;
    int i, j, k, pos;
    double basecbias;
    
    if (network != 0) {
        gtime_t basetime = timeadd(time, fmod(time2gpst(time, NULL), 30.0) * -1.0);
        int baseiod;
        if ((basebias = GetSameBiasCorrection(basetime)) == NULL) {
            return;
        }
        if ((baseiod = GetCSSRIODCBias(basebias, network)) == -1) {
            return;
        }
        if (iod != baseiod) {
            return;
        }
        trace(3, "SetCSSRBankCBias(): network=%d, tow=%.1f, iod=%d, basetow=%.1f, baseiod=%d\n",
            network, time2gpst(time, NULL), iod, time2gpst(basebias->time, NULL), baseiod);
    }

    if ((bias = GetSameBiasCorrection(time)) == NULL) {
        BiasBank[NextBias].time = time;
        BiasBank[NextBias].use = TRUE;
        bias = &BiasBank[NextBias];
        if (++NextBias >= BIASBANKNUM) {
            NextBias = 0;
        }
        memset(bias->smode, 0x00, sizeof(bias->smode));
        memset(bias->sflag, 0x00, sizeof(bias->sflag));
        memset(bias->bflag, 0x00, sizeof(bias->bflag));
        memset(bias->prn, 0x00, sizeof(bias->prn));
    }

    for (i = 0; i < MAXSAT; ++i) {
        if (network != 0) {
            if (nav->extcorr[network][i].update_cb == 1 && nav->extcorr[network][i].nsig > 0 &&
                basebias->prn[network-1][i] != 0 && (basebias->bflag[network-1] & 0x01)) {
                trace(3, "SetCSSRBankCBias(): tow=%.1f, network=%d, prn=%2d, num=%d\n", time2gpst(time, NULL),
                        network, i + 1, nav->extcorr[network][i].nsig);
                bias->prn[network-1][i] = i + 1;
                bias->bflag[network-1] |= 0x01;
                bias->udi[i] = nav->extcorr[network][i].udi[0];
                bias->iod[i] = nav->extcorr[network][i].iod[0];
                for (j = 0; j < MAXCODE; ++j) {
                    if (nav->extcorr[network][i].smode[j] != 0 && (basecbias = GetCSSRCBiasValue(basebias, network, i + 1, nav->extcorr[network][i].smode[j])) != INVALID_VALUE) {
                        if ((pos = GetAddBiasPoint(bias, network, bias->prn[network-1][i], nav->extcorr[network][i].smode[j])) != -1) {
                            if (nav->extcorr[network][i].cbias[nav->extcorr[network][i].smode[j]-1] != INVALID_VALUE) {
                                bias->cbias[network-1][i][pos] = nav->extcorr[network][i].cbias[nav->extcorr[network][i].smode[j]-1];
                                bias->cbias[network-1][i][pos] += basecbias;
                                bias->smode[network-1][i][pos] = nav->extcorr[network][i].smode[j];
                                bias->sflag[network-1][i][pos] |= 0x01;
                                trace(4, "pos=%d, mode=%2d, flag=0x%02x, pbias=%.4f\n", pos, bias->smode[network-1][i][pos],
                                        bias->sflag[network-1][i][pos], bias->cbias[network-1][i][pos]);
                            } else {
                                if (bias->sflag[network-1][i][pos] & 0x01) {
                                    bias->sflag[network-1][i][pos] ^= 0x01;
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (nav->ssr[i].update_cb == 1 && nav->ssr[i].nsig > 0) {
                trace(3, "SetCSSRBankCBias(): tow=%.1f, network=%d, prn=%2d, num=%d\n", time2gpst(time, NULL),
                        network, i + 1, nav->ssr[i].nsig);
                for (k = 0; k < CSSR_MAX_NETWORK; ++k) {
                    bias->prn[k][i] = i + 1;
                    bias->bflag[k] |= 0x01;
                    bias->udi[i] = nav->ssr[i].udi[4];
                    bias->iod[i] = nav->ssr[i].iod[4];
                    for (j = 0; j < MAXCODE; ++j) {
                        if (nav->ssr[i].smode[j] != 0) {
                            if ((pos = GetAddBiasPoint(bias, k + 1, bias->prn[k][i], nav->ssr[i].smode[j])) != -1) {
                                if (nav->ssr[i].cbias[nav->ssr[i].smode[j]-1] != INVALID_VALUE) {
                                    bias->cbias[k][i][pos] = nav->ssr[i].cbias[nav->ssr[i].smode[j]-1];
                                    bias->smode[k][i][pos] = nav->ssr[i].smode[j];
                                    bias->sflag[k][i][pos] |= 0x01;
                                    trace(4, "network=%d, pos=%d, mode=%2d, flag=0x%02x, cbias=%.4f\n", k + 1, pos,
                                            bias->smode[k][i][pos], bias->sflag[k][i][pos], bias->cbias[k][i][pos]);
                                } else {
                                    if (bias->sflag[k][i][pos] & 0x01) {
                                        bias->sflag[k][i][pos] ^= 0x01;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    trace(4, "SetCSSRBankCBias(): next=%d\n", NextBias);
}

void SetCSSRBankPBias(gtime_t time, nav_t *nav, int network)
{
    CSSRBiasBank *bias;
    int i, j, pos;

    if ((bias = GetSameBiasCorrection(time)) == NULL) {
        BiasBank[NextBias].time = time;
        BiasBank[NextBias].use = TRUE;
        bias = &BiasBank[NextBias];
        if (++NextBias >= BIASBANKNUM) {
            NextBias = 0;
        }
        memset(bias->smode, 0x00, sizeof(bias->smode));
        memset(bias->sflag, 0x00, sizeof(bias->sflag));
        memset(bias->prn, 0x00, sizeof(bias->prn));
    }

    for (i = 0; i < MAXSAT; ++i) {
        if (nav->extcorr[network][i].update_pb == 1 && nav->extcorr[network][i].nsig > 0) {
            trace(3, "SetCSSRBankPBias(): tow=%.1f, network=%d, prn=%2d, num=%d\n", time2gpst(time, NULL),
                    network, i + 1, nav->extcorr[network][i].nsig);
            bias->prn[network-1][i] = i + 1;
            bias->bflag[network-1] |= 0x02;
            bias->udi[i] = nav->extcorr[network][i].udi[1];
            bias->iod[i] = nav->extcorr[network][i].iod[1];
            for (j = 0; j < MAXCODE; ++j) {
                if (nav->extcorr[network][i].smode[j] != 0) {
                    if ((pos = GetAddBiasPoint(bias, network, bias->prn[network-1][i], nav->extcorr[network][i].smode[j])) != -1) {
                        if (nav->extcorr[network][i].pbias[nav->extcorr[network][i].smode[j]-1] != INVALID_VALUE) {
                            bias->pbias[network-1][i][pos] = nav->extcorr[network][i].pbias[nav->extcorr[network][i].smode[j]-1];
                            bias->smode[network-1][i][pos] = nav->extcorr[network][i].smode[j];
                            bias->sflag[network-1][i][pos] |= 0x02;
                            trace(4, "pos=%d, mode=%2d, flag=0x%02x, pbias=%.4f\n", pos, bias->smode[network-1][i][pos],
                                    bias->sflag[network-1][i][pos], bias->pbias[network-1][i][pos]);
                        } else {
                            if (bias->sflag[network-1][i][pos] & 0x02) {
                                bias->sflag[network-1][i][pos] ^= 0x02;
                            }
                        }
                    }
                }
            }
        }
    }
    trace(4, "SetCSSRBankPBias(): next=%d, network=%d, bflag=0x%02x\n", NextBias, network, bias->bflag[network-1]);
}

static CSSRClockBank *GetSameClockCorrection(gtime_t time)
{
    int i;

    for (i = 0; i < CLKBANKNUM; ++i) {
        if (ClockBank[i].use == TRUE && timediff(ClockBank[i].time, time) == 0.0) {
            return &ClockBank[i];
        }
    }
    return NULL;
}

static CSSRClockBank *GetCloseClockCorrection(gtime_t time, int network, double age)
{
    int search = 0;
    int pos = -1;
    int i;

    for (i = 0; i < CLKBANKNUM; ++i) {
        if (ClockBank[i].use == TRUE && (ClockBank[i].validnetwork & (1 << network))) {
            search = network;
            break;
        }
    }
    for (i = 0; i < CLKBANKNUM; ++i) {
        if (ClockBank[i].use == TRUE) {
            trace(5, "GetCloseClockCorrection(): clock=%.1f, diff=%.1f\n", time2gpst(ClockBank[i].time, NULL), timediff(ClockBank[i].time, time));
        }
        if (ClockBank[i].use == TRUE && (ClockBank[i].validnetwork & (1 << search)) && fabs(timediff(ClockBank[i].time, time)) < age) {
            if (pos != -1 && timediff(ClockBank[pos].time, ClockBank[i].time) > 0.0) {
                continue;
            }
            if (timediff(ClockBank[i].time, time) > 0.0) {
                continue;
            }
            pos = i;
        }
    }
    if (pos != -1) {
        return &ClockBank[pos];
    }
    return NULL;
}

void SetCSSRBankClock(gtime_t time, nav_t *nav, int network)
{
    CSSRClockBank *clock;
    int i;

    if ((clock = GetSameClockCorrection(time)) == NULL) {
        ClockBank[NextClock].time = time;
        ClockBank[NextClock].use = TRUE;
        clock = &ClockBank[NextClock];
        if (++NextClock >= CLKBANKNUM) {
            NextClock = 0;
        }
        memset(clock->prn, 0x00, sizeof(clock->prn));
        clock->validnetwork = 0x00;
    }
    clock->validnetwork |= (1 << network);

    for (i = 0; i < MAXSAT; ++i) {
        if (network == 0) {
            if (nav->ssr[i].update_cc == 1) {
                if (nav->ssr[i].dclk[0] != INVALID_VALUE) {
                    trace(3, "SetCSSRBankClock(): tow=%.1f, prn=%2d, c0=%.4f\n", time2gpst(time, NULL),
                            i + 1, nav->ssr[i].dclk[0]);
                    clock->prn[network][i] = i + 1;
                    clock->udi[i] = nav->ssr[i].udi[1];
                    clock->iod[i] = nav->ssr[i].iod[1];
                    clock->ssrtime[network][i] = time;
                    clock->c0[network][i] = nav->ssr[i].dclk[0];
                } else {
                    trace(3, "SetCSSRBankClock(): tow=%.1f, prn=%2d, c0=#N/A\n", time2gpst(time, NULL),
                            i + 1);
                    clock->prn[network][i] = 0;
                }
            }
        } else {
            if (nav->extcorr[network][i].update_cc == 1) {
                if (nav->extcorr[network][i].dclk != INVALID_VALUE) {
                    trace(3, "SetCSSRBankClock(): tow=%.1f, network=%d, prn=%2d, c0=%.4f\n", time2gpst(time, NULL),
                            network, i + 1, nav->extcorr[network][i].dclk);
                    clock->prn[network][i] = i + 1;
                    clock->udi[i] = nav->extcorr[network][i].udi[3];
                    clock->iod[i] = nav->extcorr[network][i].iod[3];
                    clock->ssrtime[network][i] = time;
                    clock->c0[network][i] = nav->extcorr[network][i].dclk;
                } else {
                    trace(3, "SetCSSRBankClock(): tow=%.1f, network=%d, prn=%2d, c0=#N/A\n", time2gpst(time, NULL),
                            network, i + 1);
                    clock->prn[network][i] = 0;
                }
            }
        }
    }
    trace(4, "SetCSSRBankClock(): next=%d\n", NextClock);
}

void SetCSSRLatestTrop(gtime_t time, ssrgp_t *ssrg, int network)
{
    int i, j, sat;
    
    for (i = 0; i < ssrg->ngp; ++i) {
        for (j = 0; j < ssrg->nsv[i]; ++j) {
            if (ssrg->stec[i][j] != INVALID_VALUE) {
                sat = ssrg->sat[i][j];
                LatestTrop.stec0[network-1][i][sat-1] = ssrg->stec0[i][j];
                LatestTrop.stec[network-1][i][sat-1] = ssrg->stec[i][j];
                LatestTrop.stectime[network-1][i][sat-1] = time;
                LatestTrop.prn[network-1][i][sat-1] = sat;
                trace(4, "SetCSSRLatestTrop(): network=%d, tow=%.1f, i=%2d, prn=%2d, stec=%.3f, stec0=%.3f\n", network, time2gpst(time, NULL),
                    i, sat, LatestTrop.stec[network-1][i][sat-1], LatestTrop.stec0[network-1][i][sat-1]);
            }
        }
        if (ssrg->trop_total[i] != INVALID_VALUE && ssrg->trop_wet[i] != INVALID_VALUE) {
            LatestTrop.total[network-1][i] = ssrg->trop_total[i];
            LatestTrop.wet[network-1][i] = ssrg->trop_wet[i];
            LatestTrop.troptime[network-1][i] = time;
            trace(4, "SetCSSRLatestTrop(): network=%d, tow=%.1f, i=%2d, total=%.3f, wet=%.3f\n", network, time2gpst(time, NULL),
                i, LatestTrop.total[network-1][i], LatestTrop.wet[network-1][i]);
        }
    }
    LatestTrop.gridnum[network-1] = ssrg->ngp;
}

int GetCSSRLatestTrop(double *total, double *wet, gtime_t time, int network, int index)
{
    if (index < LatestTrop.gridnum[network-1] && timediff(time, LatestTrop.troptime[network-1][index]) <= TROPVALIDAGE) {
        trace(4, "GetCSSRLatestTrop(): network=%2d, i=%2d, tow=%.1f, diff=%.1f, total=%.3f, wet=%.3f\n", network, index,
            time2gpst(time, NULL), timediff(time, LatestTrop.troptime[network-1][index]),
            LatestTrop.total[network-1][index], LatestTrop.wet[network-1][index]);
        *total = LatestTrop.total[network-1][index];
        *wet   = LatestTrop.wet[network-1][index];
        return TRUE;
    }
    trace(4, "GetCSSRLatestTrop(): nothing latest trop value, network=%2d, i=%2d, tow=%.1f\n",
        network, index, time2gpst(time, NULL));
    *total = INVALID_VALUE;
    *wet   = INVALID_VALUE;
    return FALSE;
}

int GetCSSRLatestIono(double *iono, gtime_t time, int network, int index, int sat)
{
    if (index < LatestTrop.gridnum[network-1] && LatestTrop.prn[network-1][index][sat-1] == sat &&
        timediff(time, LatestTrop.stectime[network-1][index][sat-1]) <= STECVALIDAGE) {
        if (LatestTrop.stec[network-1][index][sat-1] == INVALID_VALUE) {
            *iono = 40.3E16 / (FREQ1 * FREQ2) * LatestTrop.stec0[network-1][index][sat-1];
        } else {
            *iono = 40.3E16 / (FREQ1 * FREQ2) * LatestTrop.stec[network-1][index][sat-1];
        }
        trace(4, "GetCSSRLatestIono(): network=%2d, i=%2d, sat=%2d, tow=%.1f, diff=%.1f, iono=%.3f\n", network, index,
            sat, time2gpst(time, NULL), timediff(time, LatestTrop.stectime[network-1][index][sat-1]), *iono);
        return TRUE;
    }
    *iono = INVALID_VALUE;
    return FALSE;
}

static CSSRTropBank *GetSameTropCorrection(gtime_t time)
{
    int i;

    for (i = 0; i < TROPBANKNUM; ++i) {
        if (TropBank[i].use == TRUE && timediff(TropBank[i].time, time) == 0.0) {
            return &TropBank[i];
        }
    }
    return NULL;
}

static CSSRTropBank *GetCloseTropCorrection(gtime_t time, int network, double age)
{
    int pos = -1;
    int i;

    trace(3, "GetCloseTropCorrection(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);

    if (network >= 1 && network <= CSSR_MAX_NETWORK) {
        for (i = 0; i < TROPBANKNUM; ++i) {
            if (TropBank[i].use == TRUE) {
                trace(5, "GetCloseTropCorrection(): trop=%.1f, diff=%.1f, ngp=%d\n", time2gpst(TropBank[i].time, NULL),
                    timediff(TropBank[i].time, time), TropBank[i].gridnum[network-1]);            }
            if (TropBank[i].use == TRUE && fabs(timediff(time, TropBank[i].time)) <= age) {
                if (pos != -1 && timediff(TropBank[pos].time, TropBank[i].time) > 0.0) {
                    continue;
                }
#if 1 /* need option ? */
                if (timediff(time, TropBank[i].time) < 0.0) {
                    continue;
                }
#endif
                if (TropBank[i].gridnum[network-1] > 0) {
                    pos = i;
                }
            }
        }
        if (pos != -1) {
            trace(4, "GetCloseTropCorrection(): trop=%.1f, diff=%.1f, ngp=%d\n", time2gpst(TropBank[pos].time, NULL),
                timediff(TropBank[pos].time, time), TropBank[pos].gridnum[network-1]);
            return &TropBank[pos];
        }
    }
    return NULL;
}

void SetCSSRBankTrop(gtime_t time, ssrgp_t *ssrg, int network)
{
    CSSRTropBank *trop;
    int i, j;

    if ((trop = GetSameTropCorrection(time)) == NULL) {
        TropBank[NextTrop].time = time;
        TropBank[NextTrop].use = TRUE;
        trop = &TropBank[NextTrop];
        if (++NextTrop >= TROPBANKNUM) {
            NextTrop = 0;
        }
        memset(trop->gridnum, 0x00, sizeof(trop->gridnum));
        memset(trop->satnum, 0x00, sizeof(trop->satnum));
    }

    trop->gridnum[network-1] = ssrg->ngp;
    trace(4, "SetCSSRBankTrop(): network=%d, tow=%.1f, gridnum=%d\n", network,
        time2gpst(trop->time, NULL), trop->gridnum[network-1]);
    trace(4, "SetCSSRBankTrop(): network=%d, gridnum=%d\n", network, trop->gridnum[network-1]);
    for (i = 0; i < ssrg->ngp; ++i) {
        trop->gridpos[network-1][i][0] = ssrg->gp[i].pos[0];
        trop->gridpos[network-1][i][1] = ssrg->gp[i].pos[1];
        trop->gridpos[network-1][i][2] = ssrg->gp[i].pos[2];
        trop->total[network-1][i] = ssrg->trop_total[i];
        trop->wet[network-1][i] = ssrg->trop_wet[i];
        trop->satnum[network-1][i] = ssrg->nsv[i];
        if (trop->total[network-1][i] == INVALID_VALUE || trop->wet[network-1][i] == INVALID_VALUE) {
            GetCSSRLatestTrop(&trop->total[network-1][i], &trop->wet[network-1][i], time, network, i);
        }
        trace(4, "SetCSSRBankTrop(): i=%d, lat=%.4f, lon=%.4f, alt=%.1f, total=%.3f, wet=%.3f, satnum=%d\n", i, trop->gridpos[network-1][i][0]*R2D, trop->gridpos[network-1][i][1]*R2D,
                trop->gridpos[network-1][i][2],    trop->total[network-1][i], trop->wet[network-1][i], trop->satnum[network-1][i]);
        for (j = 0; j < trop->satnum[network-1][i]; ++j) {
            if (ssrg->stec[i][j] == INVALID_VALUE) {
                GetCSSRLatestIono(&trop->iono[network-1][i][j], time, network, i, ssrg->sat[i][j]);
            } else {
                trop->iono[network-1][i][j] = 40.3E16 / (FREQ1 * FREQ2) * ssrg->stec[i][j];
            }
            trop->prn[network-1][i][j] = ssrg->sat[i][j];
            trace(4, "SetCSSRBankTrop(): prn=%d, iono=%.3f\n", trop->prn[network-1][i][j], trop->iono[network-1][i][j]);
        }
    }
    trace(4, "SetCSSRBankTrop(): next=%d\n", NextTrop);
}

int GetCorrectionMuchTime(gtime_t time, int network)
{
    CSSROrbitBank *orbit;
    gtime_t muchtime;
    int stat = 0;

    if (!(orbit = GetCloseOrbitCorrection(time, network, 180.0))) {
        return 0;
    }
    muchtime = orbit->time;
    stat = 3;

    trace(4, "GetCorrectionMuchTime(): obstime=%.1f, stat=%d, time=%.1f\n",
            time2gpst(time, NULL), stat, time2gpst(muchtime, NULL));
    return stat;
}

void SetGridData(double *pos, int index)
{
    CurrentCSSR.stec[index].pos[0] = pos[0] * R2D;
    CurrentCSSR.stec[index].pos[1] = pos[1] * R2D;
    CurrentCSSR.stec[index].pos[2] = pos[2];
    CurrentCSSR.stec[index].n = 0;

    CurrentCSSR.zwd[index].pos[0] = (float)(pos[0] * R2D);
    CurrentCSSR.zwd[index].pos[1] = (float)(pos[1] * R2D);
    CurrentCSSR.zwd[index].pos[2] = (float)pos[2];
    CurrentCSSR.zwd[index].n = 0;
}

void InitGridIndex(CSSRBank *bank)
{
    int i;

    for (i = 0; i < RTCM_SSR_MAX_GP; ++i) {
        bank->stec[i].data = &bank->stecdata[i][0];
        bank->stec[i].nmax = MAXSAT;
        bank->zwd[i].data = &bank->zwddata[i];
        bank->zwd[i].nmax = 1;
    }
}


static int SubGetCloseCSSR(gtime_t time, int network, CSSROrbitBank **orbit, CSSRClockBank **clock, CSSRBiasBank **cbias, CSSRBiasBank **pbias, CSSRTropBank **trop, int *flag)
{
    trace(4, "SubGetCloseCSSR(): tow=%.1f, network=%d\n", time2gpst(time, NULL), network);
    *flag = FALSE;
    
    switch (GetCorrectionMuchTime(time, network)) {
    case 3:
        if (!(*orbit = GetCloseOrbitCorrection(time, network, 180.0)) || !(*trop = GetCloseTropCorrection(time, network, 30.0))) {
            return FALSE;
        }
        if (!(*cbias = GetCloseBiasCorrection2((*trop)->time, network, 30.0))) {
            if (!(*cbias = GetCloseCBiasCorrection((*orbit)->time, network, 30.0))) {
                return FALSE;
            }
            if (!(*pbias = GetClosePBiasCorrection((*orbit)->time, network, 30.0))) {
                return FALSE;
            }
        } else {
            *pbias = *cbias;
            *flag = TRUE;
        }
        break;
    default:
        return FALSE;
    }

    if (timediff((*trop)->time, (*cbias)->time) <   0.0 || timediff((*trop)->time, (*cbias)->time) >= 30.0) {
        if (*flag == TRUE) {
            trace(1, "GetCloseCSSR(): network=%d, obs=%.1f, cbias=%.1f, trop=%.1f, time is not much.\n", network,
                time2gpst(time, NULL), time2gpst((*cbias)->time, NULL), time2gpst((*trop)->time, NULL));
            return FALSE;
        }
    }
    if (timediff((*trop)->time, (*orbit)->time) < -30.0 || timediff((*trop)->time, (*orbit)->time) >= 30.0) {
        trace(1, "GetCloseCSSR(): network=%d, obs=%.1f, orbit=%.1f, trop=%.1f, time is not much.\n", network,
            time2gpst(time, NULL), time2gpst((*orbit)->time, NULL), time2gpst((*trop)->time, NULL));
        return FALSE;
    }
    if (!(*clock = GetCloseClockCorrection(time, network, 180.0))) {
        trace(1, "GetCloseCSSR(): network=%d, obs=%.1f, clock is nothing.\n", network,
            time2gpst(time, NULL));
        return FALSE;
    }
    return TRUE;
}

extern int GetCloseCSSR(gtime_t time, int network)
{
    static CSSRBank saveCSSR;
    CSSROrbitBank *orbit;
    CSSRClockBank *clock;
    CSSRBiasBank *cbias;
    CSSRBiasBank *pbias;
    CSSRTropBank *trop;
    int sis;
    int i, j;
    
    if (SubGetCloseCSSR(time, network, &orbit, &clock, &cbias, &pbias, &trop, &sis) == FALSE) {
        return FALSE;
    }
    
    trace(3, "update cssr: facility=%d, network=%d, obs=%.1f, orbit=%.1f, clock=%.1f, cbias=%.1f, pbias=%.1f, trop=%.1f, gridnum=%d\n", Facility + 1,
        network, time2gpst(time, NULL), time2gpst(orbit->time, NULL), time2gpst(clock->time, NULL), time2gpst(cbias->time, NULL),
        time2gpst(pbias->time, NULL), time2gpst(trop->time, NULL), trop->gridnum[network-1]);
    memcpy(&saveCSSR, &CurrentCSSR, sizeof(CurrentCSSR));
    memset(&CurrentCSSR, 0x00, sizeof(CurrentCSSR));
    InitGridIndex(&CurrentCSSR);
    InitGridIndex(&saveCSSR);

    for (i = 0; i < MAXSAT; ++i) {
        if (cbias->prn[network-1][i] != 0) {
            for (j = 0; j < MAXCODE; ++j) {
                if (cbias->smode[network-1][i][j] != 0 && (cbias->sflag[network-1][i][j] & 0x01) == 0x01) {
                    trace(4, "GetCloseCSSR(): cbias, network=%d, j=%d, prn=%d, mode=%d, cbias=%.6f\n", network, j,
                        cbias->prn[network-1][i], cbias->smode[network-1][i][j], cbias->cbias[network-1][i][j]);
                    CurrentCSSR.cbias[i][j] = cbias->cbias[network-1][i][j];
                    CurrentCSSR.smode[i][j] = cbias->smode[network-1][i][j];
                }
            }
        }
        if (pbias->prn[network-1][i] != 0) {
            for (j = 0; j < MAXCODE; ++j) {
                if (pbias->smode[network-1][i][j] != 0 && (pbias->sflag[network-1][i][j] & 0x02) == 0x02) {
                    trace(4, "GetCloseCSSR(): pbias, network=%d, j=%d, prn=%d, mode=%d, pbias=%.6f\n", network, j,
                        pbias->prn[network-1][i], pbias->smode[network-1][i][j], pbias->pbias[network-1][i][j]);
                    CurrentCSSR.pbias[i][j] = pbias->pbias[network-1][i][j];
                    CurrentCSSR.smode[i][j] = pbias->smode[network-1][i][j];
                }
            }
        }

        CurrentCSSR.time[i][4] = cbias->time;
        CurrentCSSR.prn[i][4] = cbias->prn[network-1][i];
        CurrentCSSR.udi[i][4] = cbias->udi[i];
        CurrentCSSR.iod[i][4] = cbias->iod[i];
        CurrentCSSR.time[i][5] = pbias->time;
        CurrentCSSR.prn[i][5] = pbias->prn[network-1][i];
        CurrentCSSR.udi[i][5] = pbias->udi[i];
        CurrentCSSR.iod[i][5] = pbias->iod[i];
        CurrentCSSR.flag[i] |= 0x04;
    }
    for (i = 0; i < MAXSAT; ++i) {
        j = (orbit->validnetwork & (1 << network) ? network: 0);
        if (orbit->prn[j][i] != 0) {
            trace(4, "GetCloseCSSR(): orbit, prn=%d, iode=%d, deph0=%f, deph1=%f, deph2=%f\n", orbit->prn[j][i],
                    orbit->iode[j][i], orbit->deph0[j][i], orbit->deph1[j][i], orbit->deph2[j][i]);
            CurrentCSSR.time[i][0] = orbit->time;
            CurrentCSSR.prn[i][0] = orbit->prn[j][i];
            CurrentCSSR.udi[i][0] = orbit->udi[i];
            CurrentCSSR.iod[i][0] = orbit->iod[i];
            CurrentCSSR.iode[i] = orbit->iode[j][i];
            CurrentCSSR.deph0[i] = orbit->deph0[j][i];
            CurrentCSSR.deph1[i] = orbit->deph1[j][i];
            CurrentCSSR.deph2[i] = orbit->deph2[j][i];
            CurrentCSSR.flag[i] |= 0x01;
        }
    }
    for (i = 0, CurrentCSSR.gridnum = trop->gridnum[network-1]; i < CurrentCSSR.gridnum; ++i) {
        for (j = 0, SetGridData(&trop->gridpos[network-1][i][0], i); j < trop->satnum[network-1][i]; ++j) {
            if (trop->iono[network-1][i][j] != INVALID_VALUE) {
                add_data_stec(&CurrentCSSR.stec[i], trop->time, trop->prn[network-1][i][j],
                        0, trop->iono[network-1][i][j], 0.0, 0.0, 0.0);
            }
        }
        add_data_trop(&CurrentCSSR.zwd[i], trop->time, trop->wet[network-1][i],
                trop->total[network-1][i], 9999, 0, 1);
    }
    for (i = 0; i < MAXSAT; ++i) {
        j = (clock->validnetwork & (1 << network) ? network: 0);
        if (clock->prn[j][i] != 0) {
            trace(4, "GetCloseCSSR(): clock, prn=%d, c0=%f\n", clock->prn[j][i], clock->c0[j][i]);
            CurrentCSSR.time[i][1] = clock->time;
            CurrentCSSR.prn[i][1] = clock->prn[j][i];
            CurrentCSSR.udi[i][1] = clock->udi[i];
            CurrentCSSR.iod[i][1] = clock->iod[i];
            CurrentCSSR.c0[i] = clock->c0[j][i];
            CurrentCSSR.flag[i] |= 0x02;
        }
    }

    CurrentCSSR.orbit_time = orbit->time;
    CurrentCSSR.clock_time = clock->time;
    CurrentCSSR.bias_time = cbias->time;
    CurrentCSSR.trop_time = trop->time;
    CurrentCSSR.facility = Facility;
    CurrentCSSR.update_time = time;
    CurrentCSSR.network = network;
    CurrentCSSR.sisadjust = sis;
    CurrentCSSR.use = TRUE;
    return CurrentCSSR.use;
}

extern int CheckGridData(gtime_t time, int network, int index)
{
    CSSROrbitBank *orbit;
    CSSRClockBank *clock;
    CSSRBiasBank *cbias;
    CSSRBiasBank *pbias;
    CSSRTropBank *trop;
    int sis;
    
    if (SubGetCloseCSSR(time, network, &orbit, &clock, &cbias, &pbias, &trop, &sis) == FALSE) {
        return FALSE;
    }
    
    trace(3, "CheckGridData(): facility=%d, network=%d, obs=%.1f, orbit=%.1f, clock=%.1f, cbias=%.1f, pbias=%.1f, trop=%.1f, gridnum=%d\n", Facility + 1,
        network, time2gpst(time, NULL), time2gpst(orbit->time, NULL), time2gpst(clock->time, NULL), time2gpst(cbias->time, NULL),
        time2gpst(pbias->time, NULL), time2gpst(trop->time, NULL), trop->gridnum[network-1]);
    
    if (index < trop->gridnum[network-1]) {
        trace(4, "grid data: network=%d, index=%d, tow=%.1f, total=%.3f, wet=%.3f, iono satnum=%d\n", network, index,
            time2gpst(trop->time, NULL), trop->total[network-1][index], trop->wet[network-1][index],
            trop->satnum[network-1][index]);
        if (trop->total[network-1][index] == INVALID_VALUE ||
            trop->wet[network-1][index] == INVALID_VALUE) {
            trace(2, "invalid grid data: network=%d, index=%d, tow=%.1f, total=%.3f, wet=%.3f\n", network, index,
                time2gpst(trop->time, NULL), trop->total[network-1][index], trop->wet[network-1][index]);
            return FALSE;
        }
        if (trop->satnum[network-1][index] < 5) {
            trace(2, "invalid grid data: network=%d, index=%d, tow=%.1f, iono satnum=%d\n", network, index,
                time2gpst(trop->time, NULL), trop->satnum[network-1][index]);
            return FALSE;
        }
        return TRUE;
    }
    return FALSE;
}

extern void UpdateGlobalCSSR(ssr_t *ssr, int sat)
{
    int i;

    if (CurrentCSSR.prn[sat-1][4] != 0) {
        for (i = ssr->nsig = 0; i < MAXCODE; ++i) {
            ssr->discontinuity[i] = 0;
            ssr->pbias[i] = 0.0;
            ssr->cbias[i] = 0.0;
            ssr->smode[i] = 0;
        }
        for (i = 0; i < MAXCODE; ++i) {
            if (CurrentCSSR.smode[sat-1][i] != 0) {
                trace(3, "UpdateGlobalCSSR(): bias correction, network=%d, tow=%.1f, sat=%d, mode=%d, cbias=%.6f, pbias=%.6f\n", CurrentCSSR.network,
                        time2gpst(CurrentCSSR.time[sat-1][4], NULL), sat, CurrentCSSR.smode[sat-1][i],
                        CurrentCSSR.cbias[sat-1][i], CurrentCSSR.pbias[sat-1][i]);
                ssr->cbias[CurrentCSSR.smode[sat-1][i]-1] = CurrentCSSR.cbias[sat-1][i];
                ssr->pbias[CurrentCSSR.smode[sat-1][i]-1] = CurrentCSSR.pbias[sat-1][i];
                ssr->smode[i] = CurrentCSSR.smode[sat-1][i];
                ++ssr->nsig;
            }
        }
        ssr->udi[4] = CurrentCSSR.udi[sat-1][4];
        ssr->udi[5] = CurrentCSSR.udi[sat-1][5];
        ssr->iod[4] = CurrentCSSR.iod[sat-1][4];
        ssr->iod[5] = CurrentCSSR.iod[sat-1][5];
        ssr->t0[4] = CurrentCSSR.time[sat-1][4];
        ssr->t0[5] = CurrentCSSR.time[sat-1][5];
        ssr->network[5] = CurrentCSSR.network;
    } else {
        memset(&ssr->t0[4], 0x00, sizeof(ssr->t0[4]));
        memset(&ssr->t0[5], 0x00, sizeof(ssr->t0[5]));
    }
    if (CurrentCSSR.prn[sat-1][0] != 0) {
        trace(3, "UpdateGlobalCSSR(): orbit correction, network=%d, tow=%.1f, sat=%d, iode=%d, deph0=%.6f\n", CurrentCSSR.network,
                time2gpst(CurrentCSSR.time[sat-1][0], NULL), sat, CurrentCSSR.iode[sat-1], CurrentCSSR.deph0[sat-1]);
        ssr->t0[0] = CurrentCSSR.time[sat-1][0];
        ssr->udi[0] = CurrentCSSR.udi[sat-1][0];
        ssr->iod[0] = CurrentCSSR.iod[sat-1][0];
        ssr->iode = CurrentCSSR.iode[sat-1];
        ssr->deph[0] = CurrentCSSR.deph0[sat-1];
        ssr->deph[1] = CurrentCSSR.deph1[sat-1];
        ssr->deph[2] = CurrentCSSR.deph2[sat-1];
        ssr->ddeph[0] = 0.0;
        ssr->ddeph[1] = 0.0;
        ssr->ddeph[2] = 0.0;
    } else {
        memset(&ssr->t0[0], 0x00, sizeof(ssr->t0[0]));
    }
    if (CurrentCSSR.prn[sat-1][1] != 0) {
        trace(3, "UpdateGlobalCSSR(): clock correction, network=%d, tow=%.1f, sat=%d, dclk=%.6f\n", CurrentCSSR.network,
                time2gpst(CurrentCSSR.time[sat-1][1], NULL), sat, CurrentCSSR.c0[sat-1]);
        ssr->t0[1] = CurrentCSSR.time[sat-1][1];
        ssr->udi[1] = CurrentCSSR.udi[sat-1][1];
        ssr->iod[1] = CurrentCSSR.iod[sat-1][1];
        ssr->dclk[0] = CurrentCSSR.c0[sat-1];
        ssr->dclk[1] = ssr->dclk[2] = 0.0;
        ssr->update = 1;
    } else {
        memset(&ssr->t0[1], 0x00, sizeof(ssr->t0[1]));
    }
}

extern void UpdateLocalCSSR(nav_t *nav)
{
    nav->stec = &CurrentCSSR.stec[0];
    nav->zwd = &CurrentCSSR.zwd[0];
}

extern int IsSISAdjust(void)
{
    return CurrentCSSR.sisadjust;
}

extern void CheckCSSRFacility(nav_t *nav, int network)
{
    static int savefacility = -1;
    static int savenetwork = -1;

    if (savefacility != -1 && savefacility != CurrentCSSR.facility) {
        trace(2, "CSSR facility changed, %d ---> %d\n", savefacility + 1, CurrentCSSR.facility + 1);
        nav->filreset = TRUE;
    }
    if (savenetwork != -1 && savenetwork != network) {
        trace(2, "CSSR network changed, %d ---> %d\n", savenetwork, network);
        if (nav->filreset != TRUE) {
            nav->ionoreset = TRUE;
        }
    }
    nav->facility = CurrentCSSR.facility + 1;
    savefacility = CurrentCSSR.facility;
    savenetwork = network;
}

extern int GetCurrentCSSRFacility(void)
{
    return(CurrentCSSR.facility + 1);
}

gtime_t GetCurrentCSSRTime(void)
{
    gtime_t time = CurrentCSSR.clock_time;

    if (timediff(CurrentCSSR.orbit_time, time) > 0.0) {
        time = CurrentCSSR.orbit_time;
    }
    if (timediff(CurrentCSSR.bias_time, time) > 0.0) {
        time = CurrentCSSR.bias_time;
    }
    return time;
}

gtime_t GetBackupCSSRTime(void)
{
    gtime_t time = BackupCSSR.clock_time;

    if (timediff(BackupCSSR.orbit_time, time) > 0.0) {
        time = BackupCSSR.orbit_time;
    }
    if (timediff(BackupCSSR.bias_time, time) > 0.0) {
        time = BackupCSSR.bias_time;
    }
    return time;
}

extern void BackupCurrentCSSR(grid_t *grid)
{
    memcpy(&BackupCSSR, &CurrentCSSR, sizeof(CSSRBank));
    memcpy(&BackupGrid, grid, sizeof(grid_t));
    InitGridIndex(&BackupCSSR);
}

extern void RestoreCurrentCSSR(grid_t *grid)
{
    memcpy(&CurrentCSSR, &BackupCSSR, sizeof(CSSRBank));
    memcpy(grid, &BackupGrid, sizeof(grid_t));
    InitGridIndex(&CurrentCSSR);
}

extern void ClearCurrentCSSR(void)
{
    memset(&BackupCSSR, 0x00, sizeof(BackupCSSR));
    InitGridIndex(&BackupCSSR);
    memset(&CurrentCSSR, 0x00, sizeof(CSSRBank));
    InitGridIndex(&CurrentCSSR);
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

    if (fp != NULL) fprintf(fp, "%d, %d, %d, %d, %d, %d", rtcm->ctype, rtcm->subtype, *tow, udi, *sync, *iod);

    if (rtcm->subtype == CSSR_TYPE_MASK) {
        cssr->iod = *iod;
        *ngnss = getbitu(rtcm->buff,i, 4); i+= 4; /* number of gnss */
        if (fp != NULL) fprintf(fp, ", %d", *ngnss);
    }

    return i;
}

/* decode mask message */
static int decode_cssr_mask(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, l, sync, tow, ngnss, iod, nsat_g=0, id, nsig, ncell=0, prn, sat[CSSR_MAX_SV];
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

    trace(2,"decode_cssr_mask: facility=%d tow=%d iod=%d\n", l6facility+1, tow, cssr->iod);

    for (k=0;k<ngnss;k++) {
        if (k != -0) {
            if (fp != NULL) fprintf(fp, ",,,,,,, ");
        } else {
            if (fp != NULL) fprintf(fp, ", ");
        }
        id = getbitu(rtcm->buff,i, 4); i+= 4; /* gnss id */
        if (fp != NULL) fprintf(fp, "%d, ", id);
        cssr->svmask[id] = (uint64_t)getbitu(rtcm->buff,i, 8)<<32; i+= 8; /* sv mask */
        if (fp != NULL) fprintf(fp, "0x%02x",  (uint32_t)(cssr->svmask[id]>>32));
        cssr->svmask[id] |= getbitu(rtcm->buff,i, 32); i+= 32; /* sv mask */
        if (fp != NULL) fprintf(fp, "%08lx, ", (uint64_t)cssr->svmask[id] & 0x00000000ffffffff);
        cssr->sigmask[id] = getbitu(rtcm->buff,i, 16); i+= 16; /* signal mask */
        if (fp != NULL) fprintf(fp, "0x%04x, ", cssr->sigmask[id]);
        cssr->cmi[id] = getbitu(rtcm->buff,i, 1); i++; /* cell mask availability */
        if (fp != NULL) fprintf(fp, "%d", cssr->cmi[id]);

        nsig = sigmask2nsig(cssr->sigmask[id]);
        nsat_g = svmask2nsatlist(cssr->svmask[id], id, sat);

        if (cssr->cmi[id]) { /* cell-mask is included */
            for (j = 0; j < nsat_g; j++) {
                cssr->cellmask[ncell] = getbitu(rtcm->buff, i, nsig); i += nsig;
                satsys(sat[j], &prn);
                if (j != 0) {
                    if (fp != NULL) fprintf(fp, ",,,,,,,,,,, %d, 0x%04x\n", prn, cssr->cellmask[ncell]);
                } else {
                    if (fp != NULL) fprintf(fp, ", %d, 0x%04x\n", prn, cssr->cellmask[ncell]);
                }
                ncell++;
            }
        } else {
            for (j = 0; j < nsat_g; j++) {
                for (l = 0; l < nsig; l++) {
                    cssr->cellmask[ncell] |= ((uint16_t)1<<(nsig-1-l));
                }
                ++ncell;
            }
            if (fp != NULL) fprintf(fp, "\n");
        }
    }
    if (ngnss <= 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }
    cssr->l6delivery = l6delivery;
    cssr->l6facility = l6facility;
    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is enough to decode the mask message */
static int check_bit_width_mask(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int k,ngnss=0,cmi=0,nsig,nsat;
    uint16_t sigmask;
    uint64_t svmask;

    i0+=20+4+1+4;
    if (i0+4>rtcm->havebit) return FALSE;
    ngnss = getbitu(rtcm->buff, i0, 4); i0+=4;

    for (k=0;k<ngnss;k++) {
        i0+=4;
        if (i0+(8+32+16+1)>rtcm->havebit) return FALSE;
        svmask = (uint64_t)getbitu(rtcm->buff, i0, 8)<<32; i0+=8;
        svmask |= (uint64_t)getbitu(rtcm->buff, i0, 32); i0+=32;
        sigmask = getbitu(rtcm->buff, i0, 16); i0+=16;
        cmi = getbitu(rtcm->buff, i0, 1); i0+=1;
        nsig = sigmask2nsig(sigmask);
        nsat = svmask2nsat(svmask);
        if (cmi) {
            if (i0+nsat*nsig>rtcm->havebit) return FALSE;
            i0+=nsat*nsig;
        }
    }
    return TRUE;
}

/* decode orbit correction */
static int decode_cssr_oc(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, iode;
    int prn, gnss;
    double udint;
    ssr_t *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    nsat = svmask2sat(cssr->svmask,sat);
    check_week_ref(rtcm, tow, ref_orbit);
    rtcm->time = gpst2time(rtcm->week_ref[ref_orbit], tow);

    trace(2,"decode_cssr_oc:   facility=%d tow=%d iod=%d\n", l6facility+1, tow, iod);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
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
        if ((gnss = sys2gnss(satsys(sat[j], &prn), NULL)) == CSSR_SYS_GAL) {
            iode = getbitu(rtcm->buff, i, 10); i+= 10; /* iode */
        } else {
            iode = getbitu(rtcm->buff, i, 8); i+= 8; /* iode */
        }

        if (j != 0) {
            if (fp != NULL) fprintf(fp, ",,,,,, %d, %d, %d, ", gnss, prn, iode);
        } else {
            if (fp != NULL) fprintf(fp, ", %d, %d, %d, ", gnss, prn, iode);
        }
        
        /* delta radial/along-track/cross-track */
        ssr->deph[0] = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;
        ssr->deph[1] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
        ssr->deph[2] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
        ssr->iode = iode;

        if (fp != NULL) {
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

        ssr->t0 [0]=rtcm->time;
        ssr->udi[0]=udint;
        ssr->iod[0]=cssr->iod;
        if (ssr->deph[0] == INVALID_VALUE || ssr->deph[1] == INVALID_VALUE || ssr->deph[2] == INVALID_VALUE) {
            ssr->deph[0] = INVALID_VALUE;
            ssr->deph[1] = INVALID_VALUE;
            ssr->deph[2] = INVALID_VALUE;
        }

        for (k=0;k<3;k++) {
            ssr->ddeph[k]=0.0;
        }

        ssr->update_oc=1;
        ssr->update=1;
        trace(4, "ssr orbit: prn=%2d, tow=%d, udi=%.1f, iod=%2d, orb=%f,%f,%f\n", sat[j], tow,
              udint, cssr->iod, ssr->deph[0],ssr->deph[1],ssr->deph[2]);
    }
    if (nsat == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }

    CheckCSSRChangedFacility(cssr->l6facility);
    SetCSSRBankOrbit(rtcm->time, &rtcm->nav, 0);
    rtcm->nav.facility = cssr->l6facility;
    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is enough to decode the orbit correction message */
static int check_bit_width_oc(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int k,sat[CSSR_MAX_SV],nsat,prn;

    i0+=37;
    if (i0>rtcm->havebit) return FALSE;
    nsat = svmask2sat(cssr->svmask, sat);
    for (k=0;k<nsat;k++) {
        i0+=(satsys(sat[k],&prn)==SYS_GAL)?51:49;
        if (i0>rtcm->havebit) return FALSE;
    }
    return TRUE;
}

/* decode clock correction */
static int decode_cssr_cc(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, iod, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, prn, gnss;
    double udint;
    ssr_t *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_clock);
    rtcm->time = gpst2time(rtcm->week_ref[ref_clock], tow);
    nsat = svmask2sat(cssr->svmask,sat);

    trace(2,"decode_cssr_cc:   facility=%d tow=%d iod=%d\n", l6facility+1, tow, iod);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
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
        gnss = sys2gnss(satsys(sat[j], &prn), NULL);
        if (j != 0) {
            if (fp != NULL) fprintf(fp, ",,,,,, %d, %d, ", gnss, prn);
        } else {
            if (fp != NULL) fprintf(fp, ", %d, %d, ", gnss, prn);
        }
        ssr->t0 [1]=rtcm->time;
        ssr->udi[1]=udint;
        ssr->iod[1]=cssr->iod;

        ssr->dclk[0] = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;

        if (ssr->dclk[0] == INVALID_VALUE) {
            trace(3, "invalid clock value: tow=%d, sat=%d, value=%d\n", tow, sat[j], ssr->dclk[0]);
            if (fp != NULL) fprintf(fp, "#N/A\n");
        } else {
            if (fp != NULL) fprintf(fp, "%f\n", (double)ssr->dclk[0]);
        }
        ssr->dclk[1] = 0.0;
        ssr->dclk[2] = 0.0;
        ssr->update_cc=1;
        ssr->update=1;
        trace(4, "ssr clock: prn=%2d, tow=%d, udi=%.1f, iod=%2d, clk=%f\n", sat[j], tow,
              udint, cssr->iod, ssr->dclk[0]);
    }
    if (nsat == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }
    CheckCSSRChangedFacility(cssr->l6facility);
    SetCSSRBankClock(rtcm->time, &rtcm->nav, 0);
    rtcm->nav.facility = cssr->l6facility;
    rtcm->nbit = i;
    return sync ? 0:10;
}

static int check_bit_width_cc(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int nsat;

    nsat = svmask2sat(cssr->svmask,NULL);
    return (i0+(12+4+1+4)+15*nsat<=rtcm->havebit);
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
    int j,k,id,*codes=NULL,sys,sys_p=-1,ofst=0,nsig_s=0,code[CSSR_MAX_SIG];
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
                case SYS_GPS: codes = (int *)codes_gps; break;
                case SYS_GLO: codes = (int *)codes_glo; break;
                case SYS_GAL: codes = (int *)codes_gal; break;
                case SYS_CMP: codes = (int *)codes_bds; break;
                case SYS_QZS: codes = (int *)codes_qzs; break;
                case SYS_SBS: codes = (int *)codes_sbs; break;
            }
            for (k=0,nsig_s=0;k<CSSR_MAX_SIG;k++) {
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

/* decode code bias message */
static int decode_cssr_cb(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, s, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, lcnt, prn, gnss;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];
    double udint;
    ssr_t *ssr=NULL;
    int first = TRUE;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_cbias);
    rtcm->time = gpst2time(rtcm->week_ref[ref_cbias], tow);
    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);

    trace(2,"decode_cssr_cb:   facility=%d tow=%d iod=%d\n", l6facility+1, tow, iod);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (k = 0; k < MAXSAT; ++k) {
        ssr = &rtcm->nav.ssr[k];
        if (rtcm->subtype == CSSR_TYPE_CB || rtcm->subtype == CSSR_TYPE_BIAS) {
            ssr->t0[4].sec = 0.0;
            ssr->t0[4].time = 0;
            ssr->udi[4]=0;
            ssr->iod[4]=0;
            ssr->update_cb=0;
            ssr->update=0;
            ssr->network[4]=0;
            ssr->satno[4]=0;
            ssr->nsig=0;
        }
        for (j = 0; j < MAXCODE; ++j) {
            if (rtcm->subtype == CSSR_TYPE_CB || rtcm->subtype == CSSR_TYPE_BIAS) {
                /* code bias */
                ssr->cbias[j] = 0.0;
                ssr->smode[j] = 0;
            }
        }
    }
    
    for (k=lcnt=0;k<nsat;k++) {
        ssr = &rtcm->nav.ssr[sat[k]-1];
        if (rtcm->subtype == CSSR_TYPE_CB || rtcm->subtype == CSSR_TYPE_BIAS) {
            ssr->t0 [4]=rtcm->time;
            ssr->udi[4]=udint;
            ssr->iod[4]=cssr->iod;
            ssr->satno[4]=sat[k];
            ssr->update_cb=1;
            ssr->update=1;
        }

        for (j=0;j<nsig[k];j++) {
            gnss = sys2gnss(satsys(sat[k], &prn), NULL);
            if (first == FALSE) {
                if (fp != NULL) fprintf(fp, ",,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
            } else {
                if (fp != NULL) fprintf(fp, ", %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                first = FALSE;
            }
            s = sig[k*CSSR_MAX_SIG+j];
            if (rtcm->subtype == CSSR_TYPE_CB || rtcm->subtype == CSSR_TYPE_BIAS) {
                /* code bias */
                ssr->cbias[s-1] = decode_sval(rtcm->buff, i, 11, 0.02); i+=11;
                if (ssr->cbias[s-1] == INVALID_VALUE) {
                    trace(3, "invalid cb value: tow=%d, sat=%d, value=%.2f\n", tow, sat[k], ssr->cbias[s-1]);
                    if (fp != NULL) fprintf(fp, "#N/A\n");
                } else {
                    if (fp != NULL) fprintf(fp, "%f\n", ssr->cbias[s-1]);
                }
                trace(4, "ssr cbias: prn=%2d, tow=%d, udi=%.1f, iod=%2d, cbias=%f\n", sat[k], tow,
                      udint, cssr->iod, ssr->cbias[s-1]);
            }
            ++lcnt;
            ssr->smode[j]=s;
        }
        ssr->nsig=nsig[k];
    }
    if (lcnt == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }
    CheckCSSRChangedFacility(cssr->l6facility);
    SetCSSRBankCBias(rtcm->time, &rtcm->nav, 0, 0);
    rtcm->nav.facility = cssr->l6facility;
    rtcm->nbit = i;
    return sync ? 0:10;
}

/* decode phase bias message */
static int decode_cssr_pb(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, gnss, prn, temp;
    int lcnt;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];
    double udint;
    int first = TRUE;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_pbias);
    rtcm->time = gpst2time(rtcm->week_ref[ref_pbias], tow);
    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);

    trace(2,"decode_cssr_pb:   facility=%d tow=%d iod=%d\n", l6facility+1, tow, iod);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    for (k=lcnt=0;k<nsat;k++) {
        gnss = sys2gnss(satsys(sat[k], &prn), NULL);

        for (j = 0; j < nsig[k]; j++) {
            if (first == FALSE) {
                if (fp != NULL) fprintf(fp, ",,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
            } else {
                if (fp != NULL) fprintf(fp, ", %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                first = FALSE;
            }
            temp = decode_sval(rtcm->buff, i, 15,0.001); i+=15;
            if (fp != NULL) {
                if (temp != INVALID_VALUE) {
                    fprintf(fp, "%f, ", (double)temp);
                } else {
                    fprintf(fp, "#N/A, ");
                }
            }
            temp = getbitu(rtcm->buff,i,2); i += 2;
            if (fp != NULL) fprintf(fp, "%d\n", temp);
            ++lcnt;
        }
    }
    if (lcnt == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }

    rtcm->nav.facility = cssr->l6facility;
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
    return i0+(12+4+1+4)+nsig_total*11<=rtcm->havebit;
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
    return i0+(12+4+1+4)+nsig_total*17<=rtcm->havebit;
}

/* code bias correction */
static int decode_cssr_bias(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, s, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat;
    int nsig[CSSR_MAX_SV], sig[CSSR_MAX_SV*CSSR_MAX_SIG], sig_p[CSSR_MAX_SV*CSSR_MAX_SIG];
    int cbflag, pbflag, netflag, network, netmask, prn, lcnt, gnss;
    double udint;
    extra_correct *ssr=NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_bias);
    rtcm->time = gpst2time(rtcm->week_ref[ref_bias], tow);

    cbflag = getbitu(rtcm->buff, i, 1); i += 1;
    if (fp != NULL) fprintf(fp, ", %d", cbflag);
    pbflag = getbitu(rtcm->buff, i, 1); i += 1;
    if (fp != NULL) fprintf(fp, ", %d", pbflag);
    netflag = getbitu(rtcm->buff, i, 1); i += 1;
    if (fp != NULL) fprintf(fp, ", %d", netflag);
    network = getbitu(rtcm->buff, i, (netflag ? 5: 0)); i += (netflag ? 5: 0);
    if (netflag == 1) {
        if (fp != NULL) fprintf(fp, ", %d", network);
    }

    nsat = svmask2sat(cssr->svmask,sat);
    sigmask2sig_p(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig_p);
    sigmask2sig(nsat, sat, cssr->sigmask, cssr->cellmask, nsig, sig);
    netmask = getbitu(rtcm->buff, i, (netflag ? nsat: 0)); i += (netflag ? nsat: 0);
    if (netflag == 1) {
        if (fp != NULL) fprintf(fp, ", 0x%08x", netmask);
    }

    trace(2,"decode_cssr_bias: facility=%d tow=%d iod=%d net=%d mask=0x%x flag=%d %d %d\n", l6facility+1, tow, iod, network, netmask, cbflag, pbflag, netflag);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }
    
    for (k = 0; k < MAXSAT; ++k) {
        ssr = &rtcm->nav.extcorr[network][k];
        if (rtcm->subtype == CSSR_TYPE_BIAS && cbflag == 1) {
            ssr->t0[0].sec = 0.0;
            ssr->t0[0].time = 0;
            ssr->udi[0]=0;
            ssr->iod[0]=0;
            ssr->update_cb=0;
            ssr->network[0]=0;
            ssr->satno[0]=0;
            ssr->facility[0]=0;
        }
        if (rtcm->subtype == CSSR_TYPE_BIAS && pbflag == 1) {
            ssr->t0[1].sec = 0.0;
            ssr->t0[1].time = 0;
            ssr->udi[1]=0;
            ssr->iod[1]=0;
            ssr->update_pb=0;
            ssr->network[1]=0;
            ssr->satno[1]=0;
            ssr->facility[1]=0;
        }
        for (j = 0; j < MAXCODE; ++j) {
            if (rtcm->subtype == CSSR_TYPE_BIAS && cbflag == 1) {
                /* code bias */
                ssr->cbias[j] = 0.0;
            }
            if (rtcm->subtype == CSSR_TYPE_BIAS && pbflag == 1) {
                /* phase bias */
                ssr->pbias[j] = 0.0;
                ssr->discontinuity[j] = 0;
            }
        }
    }
    
    for (k=lcnt=0;k<nsat;k++) {
        ssr = &rtcm->nav.extcorr[network][sat[k]-1];
        if (!((netmask>>(nsat-1-k)) & 1)) {
            continue;
        }

        if (rtcm->subtype == CSSR_TYPE_BIAS && cbflag == 1) {
            ssr->t0 [0]=rtcm->time;
            ssr->udi[0]=udint;
            ssr->iod[0]=cssr->iod;
            ssr->update_cb=1;
            ssr->network[0]=network;
            ssr->satno[0]=sat[k];
            ssr->facility[0]=cssr->l6facility+1;
        }
        if (rtcm->subtype == CSSR_TYPE_BIAS && pbflag == 1) {
            ssr->t0 [1]=rtcm->time;
            ssr->udi[1]=udint;
            ssr->iod[1]=cssr->iod;
            ssr->update_pb=1;
            ssr->network[1]=network;
            ssr->satno[1]=sat[k];
            ssr->facility[1]=cssr->l6facility+1;
        }
        gnss = sys2gnss(satsys(sat[k], &prn), NULL);
        for (j = 0; j < nsig[k]; j++) {
            if (lcnt != 0) {
                if (netflag == 1) {
                    if (fp != NULL) fprintf(fp, ",,,,,,,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                } else {
                    if (fp != NULL) fprintf(fp, ",,,,,,,,, %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
                }
            } else {
                if (fp != NULL) fprintf(fp, ", %d, %d, %d, ", gnss, prn, sig_p[k*CSSR_MAX_SIG+j]);
            }
            s = sig[k*CSSR_MAX_SIG+j];
            if (cbflag == 1) { /* code bias */
                ssr->cbias[s-1] = decode_sval(rtcm->buff, i, 11, 0.02); i+=11;
                trace(4, "ssr cbias: prn=%2d, tow=%d, udi=%.1f, iod=%2d, net=%d, s=%d, cbias=%f\n", sat[k], tow,
                      udint, cssr->iod, network, s, ssr->cbias[s-1]);
                if (fp != NULL) {
                    if (ssr->cbias[s-1] != INVALID_VALUE) {
                        fprintf(fp, "%f, ", (double)ssr->cbias[s-1]);
                    } else {
                        trace(3, "invalid cb value: tow=%d, sat=%d, value=%f\n", tow, sat[k], ssr->cbias[s-1]);
                        fprintf(fp, "#N/A, "); 
                    }
                }
            } else {
                if (fp != NULL) fprintf(fp, ", ");
            }
            if (pbflag == 1) {
                /* phase bias */
                ssr->pbias[s-1] = decode_sval(rtcm->buff, i, 15, 0.001); i+=15;
                ssr->discontinuity[s-1] = getbitu(rtcm->buff, i, 2); i+= 2;
                trace(4, "ssr pbias: prn=%2d, tow=%d, udi=%.1f, iod=%2d, net=%d, s=%d, pbias=%f\n", sat[k], tow,
                      udint, cssr->iod, network, s, ssr->pbias[s-1]);
                if (fp != NULL) { 
                    if (ssr->pbias[s-1] != INVALID_VALUE) {
                        fprintf(fp, "%f, ", (double)ssr->pbias[s-1]);
                    } else {
                        trace(3, "invalid pb value: tow=%d, sat=%d, value=%f\n", tow, sat[k], ssr->pbias[s-1]);
                        fprintf(fp, "#N/A, ");
                    }
                    fprintf(fp, "%d", ssr->discontinuity[s-1]);
                }
            } else {
                if (fp != NULL) fprintf(fp, ", ");
            }
            if (fp != NULL) fprintf(fp, "\n");
            ++lcnt;
            ssr->smode[j]=s;
        }
        ssr->nsig=nsig[k];
    }
    if (lcnt == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }
    
    if (cbflag == 1) {
        CheckCSSRChangedFacility(cssr->l6facility);
        SetCSSRBankCBias(rtcm->time, &rtcm->nav, network, cssr->iod);
    }
    if (pbflag == 1) {
        CheckCSSRChangedFacility(cssr->l6facility);
        SetCSSRBankPBias(rtcm->time, &rtcm->nav, network);
    }
    rtcm->nav.facility = cssr->l6facility;
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
    if (i0+(12+4+1+4)+3>rtcm->havebit) return FALSE;

    i0+=(12+4+1+4);
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

/* decode ura correction */
static int decode_cssr_ura(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, iod, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, gnss, prn;
    double udint;
    ssr_t *ssr = NULL;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_ura);
    rtcm->time = gpst2time(rtcm->week_ref[ref_ura], tow);
    nsat = svmask2sat(cssr->svmask,sat);

    trace(3,"decode_cssr_ura:  facility=%d tow=%d iod=%d\n", l6facility+1, tow, iod);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
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
        ssr->ura = getbitu(rtcm->buff,i, 6); i+= 6; /* ssr ura */
        ssr->update_ura=1;
        ssr->update=1;
        gnss = sys2gnss(satsys(sat[j], &prn), NULL);
        if (j != 0) {
            if (fp != NULL) fprintf(fp, ",,,,,, %d, %d, 0x%02x\n", gnss, prn, ssr->ura);
        } else {
            if (fp != NULL) fprintf(fp, ", %d, %d, 0x%02x\n", gnss, prn, ssr->ura);
        }
    }
    if (nsat == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }

    rtcm->nav.facility = cssr->l6facility;
    rtcm->nbit = i;
    return sync ? 0:10;
}

/* check if the buffer length is sufficient to decode the ura message */
static int check_bit_width_ura(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int nsat;

    nsat = svmask2sat(cssr->svmask, NULL);
    return i0+(12+4+1+4)+6*nsat<=rtcm->havebit;
}

/* decode stec correction */
static int decode_cssr_stec(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, iod, s, sync, tow, ngnss, sat[CSSR_MAX_SV], nsat, inet, a, b, gnss, prn;
    double udint;
    ssrgp_t *ssrg;
    ssrion_t *ssr_ion;
    
    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_stec);
    rtcm->time = gpst2time(rtcm->week_ref[ref_stec], tow);
    nsat = svmask2sat(cssr->svmask, sat);
    cssr->opt.stec_type = getbitu(rtcm->buff, i,  2); i+= 2; /* stec correction type */
    if (fp != NULL) fprintf(fp, ", %d, ", cssr->opt.stec_type);
    inet = getbitu(rtcm->buff, i, 5); i+= 5; /* network id */
    if (fp != NULL) fprintf(fp, "%d, ", inet);

    trace(2,"decode_cssr_stec: facility=%d tow=%d iod=%d net=%d type=%d\n", l6facility+1, tow, iod, inet, cssr->opt.stec_type);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    ssrg = &rtcm->ssrg[inet];
    ssr_ion = &rtcm->ssr_ion[inet];
    
    cssr->net_svmask[inet] = getbitu(rtcm->buff,i, nsat); i+= nsat; /* stec correction type */
    trace(4, "decode_cssr_stec: mask=0x%lx\n", cssr->net_svmask[inet]);
    if (fp != NULL) fprintf(fp, "0x%02x", (uint32_t)(cssr->net_svmask[inet]>>32));
    if (fp != NULL) fprintf(fp, "%08lx", cssr->net_svmask[inet]&0xffffffff);

    ssrg->t0 = rtcm->time;
    ssrg->udi = udint;
    ssrg->iod = iod;

    for (j=0,s=0;j<nsat;j++) {
        if ((cssr->net_svmask[inet]>>(nsat-1-j))&1) {
            gnss = sys2gnss(satsys(sat[j], &prn), NULL);
            if (s != 0) {
                if (fp != NULL) fprintf(fp, ",,,,,,,,, %d, %d, ", gnss, prn);
            } else {
                if (fp != NULL) fprintf(fp, ", %d, %d, ", gnss, prn);
            }
            ssr_ion->stec.sat[s] = sat[j];
            a = getbitu(rtcm->buff, i, 3); i+= 3;
            b = getbitu(rtcm->buff, i, 3); i+= 3;
            ssr_ion->stec.quality[s] = decode_cssr_quality_stec(a,b);
            if (fp != NULL) fprintf(fp, "0x%02x, ", (a<<3)|b);
            for (k=0;k<4;k++) ssr_ion->stec.a[s][k] = 0.0;
            
            ssr_ion->stec.a[s][0] = decode_sval(rtcm->buff, i, 14, 0.05); i+=14;
            if (fp != NULL) {
                if (ssr_ion->stec.a[s][0] != INVALID_VALUE) {
                    fprintf(fp, "%f", (double)ssr_ion->stec.a[s][0]);
                } else {
                    fprintf(fp, "#N/A");
                }
            }
            if (cssr->opt.stec_type > 0) {
                ssr_ion->stec.a[s][1] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
                if (fp != NULL) {
                    if (ssr_ion->stec.a[s][1] != INVALID_VALUE) {
                        fprintf(fp, ", %f, ", (double)ssr_ion->stec.a[s][1]);
                    } else {
                        fprintf(fp, ", #N/A, ");
                    }
                }

                ssr_ion->stec.a[s][2] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
                if (fp != NULL) {
                    if (ssr_ion->stec.a[s][2] != INVALID_VALUE) {
                        fprintf(fp, "%f", (double)ssr_ion->stec.a[s][2]);
                    } else {
                        fprintf(fp, "#N/A");
                    }
                }
            }
            if (cssr->opt.stec_type > 1) {
                ssr_ion->stec.a[s][3] = decode_sval(rtcm->buff, i, 10, 0.02); i+=10;
                if (fp != NULL) {
                    if (ssr_ion->stec.a[s][3] != INVALID_VALUE) {
                        fprintf(fp, ", %f", (double)ssr_ion->stec.a[s][3]);
                    } else {
                        fprintf(fp, ", #N/A");
                    }
                }
            }
            if (fp != NULL) fprintf(fp, "\n");
            trace(4, "decode_cssr_stec: tow=%d, sat=%d\n", tow, sat[j]);
            s++;
        }
    }
    ssr_ion->stec.network = inet;
    ssr_ion->stec.nsat = s;
    ssrg->update = 1;
    if (s == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }
    rtcm->nav.facility = cssr->l6facility;
    rtcm->nbit = i;
    trace(3, "decode_cssr_stec(): tow=%d, net=%d, bits=%d\n", tow, inet, i - i0);
    return sync ? 0:10;
}

/* check if the buffer length is sufficient to decode the stec message */
static int check_bit_width_stec(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int j,sat[CSSR_MAX_SV],nsat,stec_type,slen=0,nsat_local=0;
    uint64_t net_svmask;
    const int slen_t[4] = {20,44,54,0};

    nsat = svmask2sat(cssr->svmask, sat);
    if (i0+(12+4+1+4)+2+5+nsat>rtcm->havebit) return FALSE;

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

/* decode grid correction */
static int decode_cssr_grid(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, k, ii, s, sync, iod, tow, ngnss, sat[CSSR_MAX_SV], nsat, sz;
    int trop_type, sz_idx, inet, a,b, hs, wet, gnss, prn;
    double udint, stec0, dlat, dlon, dstec;
    nav_t *nav = &rtcm->nav;
    ssrgp_t *ssrg;
    ssrion_t *ssr_ion;
    int valid;

    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_grid);
    rtcm->time = gpst2time(rtcm->week_ref[ref_grid], tow);
    nsat = svmask2sat(cssr->svmask,sat);

    trop_type = getbitu(rtcm->buff, i, 2); i+=2;  /* troposphere correction type */
    if (fp != NULL) fprintf(fp, ", %d, ", trop_type);
    sz_idx = getbitu(rtcm->buff, i, 1); i++; /* stec range */
    if (fp != NULL) fprintf(fp, "%d, ", sz_idx);
    inet = getbitu(rtcm->buff, i, 5); i+=5; /* network id */
    if (fp != NULL) fprintf(fp, "%d, ", inet);

    ssrg = &rtcm->ssrg[inet];
    ssr_ion = &rtcm->ssr_ion[inet];

    cssr->net_svmask[inet] = getbitu(rtcm->buff, i, nsat); i+= nsat; /* stec correction type */
    if (fp != NULL) fprintf(fp, "0x%02x", (uint32_t)(cssr->net_svmask[inet]>>32));
    if (fp != NULL) fprintf(fp, "%08lx, ", cssr->net_svmask[inet]&0xffffffff);
    a = getbitu(rtcm->buff, i, 3); i+= 3;
    b = getbitu(rtcm->buff, i, 3); i+= 3;
    if (fp != NULL) fprintf(fp, "0x%02x, ", (a<<3)|b);
    ssrg->ngp = getbitu(rtcm->buff,i, 6); i+= 6;
    if (fp != NULL) fprintf(fp, "%d", ssrg->ngp);
    ssrg->t0 = rtcm->time;
    ssrg->quality = decode_cssr_quality_trop(a, b);
    ssrg->network = inet;

    trace(2,"decode_cssr_grid: facility=%d tow=%d iod=%d net=%d trop=%d sz_idx=%d ngp=%d\n", l6facility+1, tow, iod, inet, trop_type, sz_idx, ssrg->ngp);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
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
        if (j != 0) {
            if (fp != NULL) fprintf(fp, ",,,,,,,,,,,, %d, ", j+1);
        } else {
            if (fp != NULL) fprintf(fp, ", %d, ", j+1);
        }
        valid=1;
        switch (trop_type) {
            case 0: break;
            case 1:
                hs = getbits(rtcm->buff, i, 9); i+=9;
                wet = getbits(rtcm->buff, i, 8); i+= 8;
                if (hs==(-P2_S9_MAX-1)) {
                    trace(2, "trop(hs) is invalid: tow=%d, inet=%d, grid=%d, hs=%d\n",
                        tow, inet, j, hs);
                    valid=0;
                    if (fp != NULL) fprintf(fp, "#N/A, ");
                } else {
                    if (fp != NULL) fprintf(fp, "%f, ", (double)hs*0.004);
                }
                if (wet==(-P2_S8_MAX-1)) {
                    trace(2, "trop(wet) is invalid: tow=%d, inet=%d, grid=%d, wet=%d\n",
                        tow, inet, j, wet);
                    valid=0;
                    if (fp != NULL) fprintf(fp, "#N/A");
                } else {
                    if (fp != NULL) fprintf(fp, "%f", (double)wet*0.004);
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
                gnss = sys2gnss(satsys(sat[k], &prn), NULL);
                if (ii != 0) {
                    if (fp != NULL) fprintf(fp, ",,,,,,,,,,,,,,, %d, %d, ", gnss, prn);
                } else {
                    if (fp != NULL) fprintf(fp, ", %d, %d, ", gnss, prn);
                }
                dstec = decode_sval(rtcm->buff, i, sz, 0.04); i+=sz;
                stec0 = ssr_ion->stec.a[ii][0] + ssr_ion->stec.a[ii][1]*dlat +
                        ssr_ion->stec.a[ii][2]*dlon +
                        ssr_ion->stec.a[ii][3]*dlat*dlon;
                if (dstec == INVALID_VALUE) {
                    trace(2, "dstec is invalid: tow=%d, inet=%d, grid=%d, sat=%d, dstec=%d\n",
                        tow, inet, j, sat[k], dstec);
                    ssrg->stec[j][s] = INVALID_VALUE;
                    if (fp != NULL) fprintf(fp, "#N/A\n");
                } else {
                    if (fp != NULL) fprintf(fp, "%f\n", (double)dstec);
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
        if (ii == 0) {
            if (fp != NULL) fprintf(fp, "\n");
        }
    }
    if (ssrg->ngp == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }

    ssrg->update = 1;
    nav->updateac=1;

    CheckCSSRChangedFacility(cssr->l6facility);
    SetCSSRLatestTrop(ssrg->t0, ssrg, inet);
    SetCSSRBankTrop(ssrg->t0, ssrg, inet);
    rtcm->nav.facility = cssr->l6facility;
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

/* decode orbit/clock combination message */
static int decode_cssr_combo(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp)
{
    int i, j, sync, iod, tow, ngnss, sat[CSSR_MAX_SV], nsat, iode;
    int flg_orbit, flg_clock, flg_net, net_svmask, netid, s, gnss, prn;
    double udint;
    extra_correct *ssr=NULL;
    
    i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp);
    check_week_ref(rtcm, tow, ref_combined);
    rtcm->time = gpst2time(rtcm->week_ref[ref_combined], tow);
    nsat = svmask2sat(cssr->svmask, sat);
    
    trace(2, "decode_cssr_combo:facility=%d tow=%d iod=%d\n", l6facility+1, tow, iod);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }
    
    flg_orbit = getbitu(rtcm->buff, i, 1); i += 1;
    if (fp != NULL) fprintf(fp, ", %d", flg_orbit);
    flg_clock = getbitu(rtcm->buff, i, 1); i += 1;
    if (fp != NULL) fprintf(fp, ", %d", flg_clock);
    flg_net = getbitu(rtcm->buff, i, 1); i += 1;
    if (fp != NULL) fprintf(fp, ", %d", flg_net);
    netid = getbitu(rtcm->buff, i, 5); i += 5;
    if (fp != NULL) fprintf(fp, ", %d", netid);
    net_svmask = getbitu(rtcm->buff, i, nsat); i += nsat;
    if (fp != NULL) fprintf(fp, ", 0x%04x", net_svmask & 0xffff);
    
    for (j = 0; j < MAXSAT; ++j) {
        ssr = &rtcm->nav.extcorr[netid][j];
        if (flg_orbit == 1) {
            ssr->t0[2].sec = 0.0;
            ssr->t0[2].time = 0;
            ssr->udi[2] = 0;
            ssr->iod[2] = 0;
            ssr->network[2] = 0;
            ssr->satno[2] = 0;
            ssr->facility[2] = 0;
            ssr->update_oc = 0;
            ssr->iode = 0;
            ssr->deph[0] = 0.0;
            ssr->deph[1] = 0.0;
            ssr->deph[2] = 0.0;
        }
        if (flg_clock == 1) {
            ssr->t0[3].sec = 0.0;
            ssr->t0[3].time = 0;
            ssr->udi[3] = 0;
            ssr->iod[3] = 0;
            ssr->network[3] = 0;
            ssr->satno[3] = 0;
            ssr->facility[3] = 0;
            ssr->update_cc = 0;
            ssr->dclk = 0.0;
        }
    }
    
    for (j = s = 0; j < nsat; ++j) {
        if ((net_svmask >> (nsat - 1 - j)) & 1) {
            ssr = &rtcm->nav.extcorr[netid][sat[j]-1];
            gnss = sys2gnss(satsys(sat[j], &prn), NULL);
            if (s != 0) {
                if (fp != NULL) fprintf(fp, ",,,,,,,,,,, %d, %d", gnss, prn);
            } else {
                if (fp != NULL) fprintf(fp, ", %d, %d", gnss, prn);
            }

            if (flg_orbit == 1) {
                if (satsys(sat[j], NULL) == SYS_GAL) {
                    iode = getbitu(rtcm->buff, i, 10); i += 10; /* iode */
                    if (fp != NULL) fprintf(fp, ", %d", iode);
                } else {
                    iode = getbitu(rtcm->buff, i,  8); i +=  8; /* iode */
                    if (fp != NULL) fprintf(fp, ", %d", iode);
                }
                /* delta radial,along-track,cross-track */
                ssr->deph[0] = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;
                ssr->deph[1] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
                ssr->deph[2] = decode_sval(rtcm->buff, i, 13, 0.0064); i+=13;
                if (fp != NULL) {
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
                }
                ssr->t0 [2] = rtcm->time;
                ssr->udi[2] = udint;
                ssr->iod[2] = cssr->iod;
                ssr->network[2] = netid;
                ssr->satno[2] = sat[j];
                ssr->facility[2] = cssr->l6facility + 1;
                ssr->update_oc = 1;
                
                ssr->iode = iode;
                if (ssr->deph[0] == INVALID_VALUE || ssr->deph[1] == INVALID_VALUE || ssr->deph[2] == INVALID_VALUE) {
                    trace(3, "invalid orbit value: tow=%d, sat=%d, value=%f %f %f\n", tow, sat[j], ssr->deph[0], ssr->deph[1], ssr->deph[2]);
                    ssr->deph[0] = INVALID_VALUE;
                    ssr->deph[1] = INVALID_VALUE;
                    ssr->deph[2] = INVALID_VALUE;
                }
                trace(4, "combined orbit: network=%d, tow=%d, sat=%d, iode=%d, deph=%f, %f, %f\n", netid, tow, sat[j], ssr->iode,
                      ssr->deph[0], ssr->deph[1], ssr->deph[2]);
            } else {
                if (fp != NULL) fprintf(fp, ",,,,");
            }
            if (flg_clock == 1) {
                ssr->dclk = decode_sval(rtcm->buff, i, 15, 0.0016); i+=15;
                
                ssr->t0 [3] = rtcm->time;
                ssr->udi[3] = udint;
                ssr->iod[3] = cssr->iod;
                ssr->network[3] = netid;
                ssr->satno[3] = sat[j];
                ssr->facility[3] = cssr->l6facility + 1;
                ssr->update_cc = 1;
                
                if (ssr->dclk == INVALID_VALUE) {
                    trace(3, "invalid clock value: tow=%d, sat=%d, value=%d\n", tow, sat[j], ssr->dclk);
                    if (fp != NULL) fprintf(fp, ", #N/A");
                } else {
                    if (fp != NULL) fprintf(fp, ", %f", (double)ssr->dclk);
                }
                trace(4, "combined clock: network=%d, tow=%d, sat=%d, dclk=%f\n", netid, tow, sat[j], ssr->dclk);
            }
            if (fp != NULL) fprintf(fp, "\n");
            ++s;
        }
    }
    if (s == 0) {
        if (fp != NULL) fprintf(fp, "\n");
    }
    
    rtcm->nav.separation |= (1 << (netid - 1));
    CheckCSSRChangedFacility(cssr->l6facility);

    if (flg_orbit == 1) {
        SetCSSRBankOrbit(rtcm->time, &rtcm->nav, netid);
    }
    if (flg_clock == 1) {
        SetCSSRBankClock(rtcm->time, &rtcm->nav, netid);
    }
    rtcm->nav.facility = cssr->l6facility;
    rtcm->nbit = i;
    return sync ? 0: 10;
}

static int check_bit_width_combo(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int sat[CSSR_MAX_SV],nsat,j,flg_orbit,flg_clock,flg_net,sz;
    uint64_t net_svmask=0;

    nsat = svmask2sat(cssr->svmask,sat);
    if (i0+(12+4+1+4)+3>rtcm->havebit) return FALSE;

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
        if ((net_svmask>>(nsat-1-j))&1) {
            if (flg_orbit) {
                sz = (satsys(sat[j],NULL)==SYS_GAL) ? 10:8;
                i0+=sz+41;
                if (i0>rtcm->havebit) return FALSE;
            }
            if (flg_clock) {
                if ((i0+=15)>rtcm->havebit) return FALSE;
            }
        }
    }
    return TRUE;
}

static int decode_cssr_atmos(rtcm_t *rtcm, cssr_t *cssr, int i0, int header, FILE *fp1, FILE *fp2)
{
    int i, j, k, s, sync, tow, iod, ngnss, gnss, sat[CSSR_MAX_SV], nsat, sz_idx, sz, prn;
    int trop_flag, stec_flag, trop_type=-1, stec_type=-1, inet, a, b, quality;
    double total[CSSR_MAX_GP], wet[CSSR_MAX_GP], stec[CSSR_MAX_GP][CSSR_MAX_SV];
    double udint, stec0, ct[6]={0}, ci[6]={0}, trop_ofst, trop_residual, dlat, dlon, dstec;
    nav_t *nav = &rtcm->nav;
    ssrgp_t *ssrg;
    const double dstec_lsb_t[4] = {0.04,0.12,0.16,0.24};
    const int dstec_sz_t[4] = {4,4,5,7};
    const int ct_num[3]={1, 3, 4}, ci_num[4]={1, 3, 4, 6};

    if ((i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp1)) == -1) {
        if (fp1 != NULL) return -1;
    }
    if ((i = decode_cssr_head(rtcm, cssr, &sync, &tow, &iod, NULL, &udint, &ngnss, i0, header, fp2)) == -1) {
        if (fp2 != NULL) return -1;
    }
    check_week_ref(rtcm, tow, ref_atmospheric);
    rtcm->time = gpst2time(rtcm->week_ref[ref_atmospheric], tow);
    nsat = svmask2sat(cssr->svmask, sat);
    
    trop_flag = getbitu(rtcm->buff, i, 2); i += 2;  /* troposphere correction availability */
    if (fp1 != NULL) fprintf(fp1, ", %d", trop_flag);
    stec_flag = getbitu(rtcm->buff, i, 2); i += 2;  /* stec correction availability */
    if (fp1 != NULL) fprintf(fp1, ", %d", stec_flag);
    inet = getbitu(rtcm->buff, i, 5); i += 5;       /* network id */
    if (fp1 != NULL) fprintf(fp1, ", %d", inet);
    if (fp2 != NULL) fprintf(fp2, ", %d", inet);
    
    trace(2, "decode_cssr_atmos:facility=%d tow=%d iod=%d net=%d\n", l6facility+1, tow, iod, inet);
    if (cssr->l6facility != l6facility) {
        trace(2, "cssr: facility mismatch: tow=%d mask_facility=%d subtype=%d facility=%d\n", tow, cssr->l6facility, rtcm->subtype, l6facility);
        return -1;
    }
    if (cssr->iod != iod) {
        trace(2, "cssr: iod mismatch: tow=%d mask_iod=%d subtype=%d iod=%d\n", tow, cssr->iod, rtcm->subtype, iod);
        return -1;
    }

    ssrg = &rtcm->ssrg[inet];
    
    ssrg->ngp = getbitu(rtcm->buff, i, 6); i += 6;
    if (fp1 != NULL) fprintf(fp1, ", %d", ssrg->ngp);
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

	if (trop_flag != 0) {
        a = getbitu(rtcm->buff, i, 3); i += 3;
        b = getbitu(rtcm->buff, i, 3); i += 3;
        ssrg->quality = decode_cssr_quality_trop(a, b);
        if (fp1 != NULL) fprintf(fp1, ", 0x%02x", (a<<3)|b);
	}
    
    if ((trop_flag&0x01) == 0x01) {
        trop_type = getbitu(rtcm->buff, i, 2); i += 2;
        if (fp1 != NULL) fprintf(fp1, ", %d", trop_type);
        for (k=0;k<4;k++) ct[k] = 0.0;
        ct[0] = decode_sval(rtcm->buff, i, 9, 0.004); i+=9;
        if (trop_type>0) {
            ct[1] = decode_sval(rtcm->buff, i, 7, 0.002); i+=7;
            ct[2] = decode_sval(rtcm->buff, i, 7, 0.002); i+=7;
        }
        if (trop_type>1) {
            ct[3] = decode_sval(rtcm->buff, i, 7, 0.001); i+=7;
        }
        for (k=0; k<4; k++) {
            if (fp1 == NULL) continue;
            if (k < ct_num[trop_type]) {
                fprintf(fp1, ", %.3f", ct[k]);
            } else {
                fprintf(fp1, ", ");
            }
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

    if ((trop_flag&0x02) == 0x02) {
        sz_idx = getbitu(rtcm->buff, i, 1); i += 1;
        if (fp1 != NULL) fprintf(fp1, ", %d", sz_idx);
        trop_ofst = getbitu(rtcm->buff, i, 4) * 0.02; i += 4;
        if (fp1 != NULL) fprintf(fp1, ", %.2f", trop_ofst);
        trace(3, "decode_cssr_atmos: network=%d, tow=%d, trop=0x%02x, stec=0x%02x, ngp=%d, trop_type=%d, ct=%.3f %.3f %.3f %.3f, sz_idx=%d, offset=%.3f\n",
            ssrg->network, tow, trop_flag, stec_flag, ssrg->ngp, trop_type, ct[0], ct[1], ct[2], ct[3], sz_idx, trop_ofst);
        sz = (sz_idx==0) ? 6:8;
        
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
    
    if (stec_flag) {
        cssr->net_svmask[inet] = getbitu(rtcm->buff, i, nsat); i += nsat; /* stec correction type */
        if (fp1 != NULL) fprintf(fp1, ", 0x%lx", cssr->net_svmask[inet]);
        trace(4, "decode_cssr_atmos: mask=0x%lx\n", cssr->net_svmask[inet]);
        
        for (j = s = 0; j < nsat; ++j) {
            if (!((cssr->net_svmask[inet] >> (nsat - 1 - j)) & 1)) {
                continue;
            }
            gnss = sys2gnss(satsys(sat[j], &prn), NULL);
            if (s != 0) {
                if (fp1 != NULL) fprintf(fp1, "\n,,,,,,,,,,,,,,,,,,, %d, %d", gnss, prn);
            } else {
                if (fp1 != NULL) fprintf(fp1, ", %d, %d", gnss, prn);
            }
            a = getbitu(rtcm->buff, i, 3); i+=3;
            b = getbitu(rtcm->buff, i, 3); i+=3;
            quality = decode_cssr_quality_stec(a,b);
            if (fp1 != NULL) fprintf(fp1, ", 0x%02x", (a<<3)|b);
            stec_type = getbitu(rtcm->buff, i, 2); i += 2;
            if (fp1 != NULL) fprintf(fp1, ", %d", stec_type);
            
            for (k=0;k<6;k++) ci[k]=0.0;
            ci[0] = decode_sval(rtcm->buff, i, 14, 0.05); i+=14;
            if (stec_type>0) {
                ci[1] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
                ci[2] = decode_sval(rtcm->buff, i, 12, 0.02); i+=12;
            }
            if (stec_type>1) {
                ci[3] = decode_sval(rtcm->buff, i, 10, 0.02); i+=10;
            }
            if (stec_type>2) {
                ci[4] = decode_sval(rtcm->buff, i, 8, 0.005); i+=8;
                ci[5] = decode_sval(rtcm->buff, i, 8, 0.005); i+=8;
            }            
            for (k=0; k<6; k++) {
                if (fp1 == NULL) continue;
                if (k < ci_num[stec_type]) {
                    if (k > 3) {
                        fprintf(fp1, ", %.3f", ci[k]);
                    } else {
                        fprintf(fp1, ", %.2f", ci[k]);
                    }
                } else {
                    fprintf(fp1, ", ");
                }
            }
            
            sz_idx = getbitu(rtcm->buff,i,2); i+=2;
            if (fp1 != NULL) fprintf(fp1, ", %d", sz_idx);
            trace(3, "decode_cssr_atmos: stec_type=%d, ct=%.2f %.2f %.2f %.2f %.2f %.2f, sz_idx=%d\n",
                stec_type, ci[0], ci[1], ci[2], ci[3], ci[4], ci[5], sz_idx);
            
            for (k = 0; k < ssrg->ngp; ++k) {
                dlat = (ssrg->gp[k].pos[0] - ssrg->gp[0].pos[0]) * R2D;
                dlon = (ssrg->gp[k].pos[1] - ssrg->gp[0].pos[1]) * R2D;
                
                dstec = decode_sval(rtcm->buff, i, dstec_sz_t[sz_idx], dstec_lsb_t[sz_idx]);
                i += dstec_sz_t[sz_idx];
                trace(5,"sz_idx=%d, dstec=%f\n",dstec_sz_t[sz_idx],dstec);
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
        
        for (k = 0; k < ssrg->ngp; ++k) {
            trace(4, "decode_cssr_atmos: grid=%d, nsv=%d\n", k + 1, ssrg->nsv[k]);
        }
    }

    if (fp2 != NULL) fprintf(fp2, ", 0x%lx", cssr->net_svmask[inet]);
    if (fp2 != NULL) fprintf(fp2, ", %d", ssrg->ngp);

    for (j = 0; j < ssrg->ngp; ++j) {
        if (j != 0) {
            if (fp2 != NULL) fprintf(fp2, "\n,,,,,,,,, %d", j + 1);
        } else {
            if (fp2 != NULL) fprintf(fp2, ", %d", j + 1);
        }

        if (trop_flag != 0) {
            if (total[j] != INVALID_VALUE) {
                if (fp2 != NULL) fprintf(fp2, ", %.3f", total[j]);
            } else {
                if (fp2 != NULL) fprintf(fp2, ", #N/A");
            }
            if (wet[j] != INVALID_VALUE) {
                if (fp2 != NULL) fprintf(fp2, ", %.3f", wet[j]);
            } else {
                if (fp2 != NULL) fprintf(fp2, ", #N/A");
            }
        } else {
            if (fp2 != NULL) fprintf(fp2, ", ,");
        }

        if (stec_flag != 0) {
            for (k = s = 0; k < nsat; ++k) {
                if (!((cssr->net_svmask[inet]  >> (nsat - 1 - k)) & 1)) {
                    continue;
                }

                gnss = sys2gnss(satsys(sat[k], &prn), NULL);
                if (s != 0) {
                    if (fp2 != NULL) fprintf(fp2, "\n,,,,,,,,,,,, %d, %d", gnss, prn);
                } else {
                    if (fp2 != NULL) fprintf(fp2, ", %d, %d", gnss, prn);
                }

                if (stec[j][s] != INVALID_VALUE) {
                    if (fp2 != NULL) fprintf(fp2, ", %.4f", stec[j][s]);
                } else {
                    if (fp2 != NULL) fprintf(fp2, ", #N/A");
                }
                ++s;
            }
        }
    }

    if (fp1 != NULL) fprintf(fp1, "\n");
    if (fp2 != NULL) fprintf(fp2, "\n");
    
    ssrg->update = 1;
    nav->updateac=1;
    
    CheckCSSRChangedFacility(cssr->l6facility);
    SetCSSRLatestTrop(ssrg->t0, ssrg, inet);
    SetCSSRBankTrop(ssrg->t0, ssrg, inet);
    
    rtcm->nav.facility = cssr->l6facility;
    rtcm->nbit = i;
    trace(3, "decode_cssr_atmos(): tow=%d, net=%d, bits=%d\n", tow, inet, i - i0);
    return sync ? 0: 10;
}

static int check_bit_width_atmos(rtcm_t *rtcm, cssr_t *cssr, int i0)
{
    int flg_trop,flg_stec,trop_type,stec_type,ngp,sz_idx,sz,j,nsat;
    uint64_t net_svmask;
    const int dstec_sz_t[4] = {4,4,5,7};
    const int trop_sz_t[3] = {9,23,30};
    const int stec_sz_t[4] = {14,38,48,64};
    
    nsat = svmask2sat(cssr->svmask,NULL);
    if (i0+(12+4+1+4)+2+2+5+6>rtcm->havebit) return FALSE;
    i0+=21;

    flg_trop = getbitu(rtcm->buff, i0, 2); i0+=2;
    flg_stec = getbitu(rtcm->buff, i0, 2); i0+=2;
    i0 += 5;
    ngp = getbitu(rtcm->buff, i0, 6); i0+=6;
    
    if (flg_trop) {
        if (i0+8>rtcm->havebit) return FALSE;
        i0+=6;
        trop_type = getbitu(rtcm->buff, i0, 2); i0+=2;
        sz = trop_sz_t[trop_type];
        if (i0+sz+5>rtcm->havebit) return FALSE;
        i0+=sz;
        sz_idx = getbitu(rtcm->buff, i0, 1); i0+=1;
        i0+=4;
        sz = (sz_idx==0)?6:8;
        if (i0+sz*ngp>rtcm->havebit) return FALSE;
        i0+=sz*ngp;
    }
    
    if (flg_stec) {
        if (i0+nsat>rtcm->havebit) return FALSE;
        net_svmask = getbitu(rtcm->buff, i0, nsat); i0+=nsat;
        for (j=0;j<nsat;j++) {
            if (!((net_svmask>>(nsat-1-j))&1)) continue;
            if (i0+8>rtcm->havebit) return FALSE;
            i0+=6;
            stec_type = getbitu(rtcm->buff, i0, 2); i0+=2;
            sz = stec_sz_t[stec_type];
            if (i0+sz+2>rtcm->havebit) return FALSE;
            i0+=sz;
            sz_idx = getbitu(rtcm->buff, i0, 2); i0+=2;
            i0+=ngp*dstec_sz_t[sz_idx];
            if (i0>rtcm->havebit) return FALSE;
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
        fgets(buff, sizeof(buff), fp);
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
    static int startdecode = FALSE;
    static cssr_t _cssr = {0,};
    static int savefacility = -1;
    static int savedelivery = -1;

    cssr_t *cssr = &_cssr;
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
    if (rtcm->subtype != 1 && startdecode == FALSE) {
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
    case 0:
        trace(1, "invalid process: frame=%02x, havebit=%d, nbit=%d\n", *frame, rtcm->havebit, rtcm->nbit);
        return 0;
    }
    
    trace(4, "cssr: frame=0x%02x, ctype=%d, subtype=%d\n", *frame, rtcm->ctype, rtcm->subtype);
    tow = time2gpst(timeget(), &week);
    if (savefacility != l6facility) {
        if (savefacility != -1) {
            trace(1, "L6 data: change facility, week=%d, tow=%.2f, %d(%d) ---> %d(%d)\n", week, tow,
                  (savefacility == -1 ? 0: savefacility+1), (savefacility == -1 ? 0: savedelivery),
                  l6facility+1, l6delivery);
        } else {
            trace(1, "L6 data: change facility, week=%d, tow=%.2f,        ---> %d(%d)\n", week, tow,
                  l6facility+1, l6delivery);
        }
        savedelivery = l6delivery;
        savefacility = l6facility;
    }
    
    switch (rtcm->subtype) {
        case CSSR_TYPE_MASK:
            ret=decode_cssr_mask(rtcm, cssr, i, head, ofp[0]);
            if (startdecode == FALSE) {
                trace(1, "start CSSR decoding: week=%d, tow=%.1f\n", week, tow);
                startdecode = TRUE;
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

/* decode cssr messages in the QZS L6 subframe */
extern int input_cssr(rtcm_t *cssr, unsigned char data, uint8_t *frame)
{
    static uint32_t preamble = 0;
    static uint64_t data_p = 0;
    uint8_t prn, msgid, alert;
    static int nframe = 0;
    static unsigned char buff[BLEN_MSG];
    static int decode_start = 0;

    trace(5,"input_cssr: data=%02x\n",data);

    /* synchronize frame */
    if (cssr->nbyte==0) {
        preamble = (preamble << 8) | data;
        data_p = (data_p << 8) | data;
        if (preamble != L6FRMPREAMB) {
            return 0;
        }
        preamble = 0;
        buff[cssr->nbyte++]=(L6FRMPREAMB>>24) & 0xff;
        buff[cssr->nbyte++]=(L6FRMPREAMB>>16) & 0xff;
        buff[cssr->nbyte++]=(L6FRMPREAMB>>8) & 0xff;
        buff[cssr->nbyte++]=data;
        return 0;
    }
    buff[cssr->nbyte++]=data;
    cssr->len = BLEN_MSG;

    if (cssr->nbyte<cssr->len) return 0;
    cssr->nbyte=0;

    prn = buff[4];
    msgid = buff[5];
    alert = (buff[6]>>7) & 0x1;

    l6delivery = prn;
    l6facility = (msgid & 0x18) >> 3;

    if (msgid & 0x01) { /* Subframe indicator */
        if (decode_start == 0) {
            trace(1, "CSSR frame first recieve: tow=%.1f\n", time2gpst(timeget(), NULL));
        }
        cssr->havebit = 0;
        decode_start = 1;
        cssr->nbit = 0;
        *frame = 0;
        nframe = 0;
    } else if (nframe >= 5) {
        return 0;
    }
    if (decode_start == 1) {
        int i = 1695 * nframe, j;

        setbitu(cssr->buff, i, 7, buff[6] & 0x7f); i+=7;
        for (j = 0; j < 211; j++) {
            setbitu(cssr->buff, i, 8, buff[7+j]); i+=8;
        }
        cssr->havebit += 1695;
        *frame |= (1<<nframe);
        nframe++;
    }
    return 0;
}


/* decode cssr messages from file stream ---------------------------------------------*/
extern int input_cssrf(rtcm_t *cssr, FILE *fp, FILE **ofp)
{
    static uint8_t frame = 0;
    int i,data=0,ret;

    trace(4,"input_cssrf: data=%02x\n",data);

    for (i=0;i<4096;i++) {
        if ((ret=decode_qzs_msg(cssr, 0, &frame, ofp))) return ret;
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_cssr(cssr, (unsigned char)data, &frame))) return ret;
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
        return -1;
    }

    while (1) {
        if ((ret = input_cssrf(&rtcm, fp, ofp)) < -1) {
            break;
        }
    }
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

