/*------------------------------------------------------------------------------
*  ssr2obs : convert Compact SSR to RINEX3 or RTCM3 MSM
*
*  Copyright (C) 2007- by T.TAKASU, All rights reserved.
*  Copyright (C) 2015- by Mitsubishi Electric Corp., All rights reserved.
*
* history : 2019/02/01  1.0  new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "cssr.h"

/* constants and macros ------------------------------------------------------*/

#define PROGNAME    "ssr2obs"           /* program name */
#define PROG_VER    "1.0"             /* program version */

#define OSR_RCORR   1                 /* output mode: range corrections */
#define OSR_RTCM3   2                 /* output mode: RTCM 3 MSM */
#define OSR_RINEX   3                 /* output mode: RINEX OBS */
#define OSR_COSR    4                 /* output mode: correction in OSR */
#define OSR_L6DUMP  5                 /* output mode: dump L6 messages */

#define OSR_SYS     (SYS_GPS|SYS_QZS|SYS_GAL) /* navigation systems */
#define OSR_NFREQ   4                 /* number of frequencies */
#define OSR_ELMASK  0.0               /* elevation mask (deg) */
#define OSR_RNXVER  3.02              /* RINEX version */
#define OSR_MARKER  "CLAS_FOR_VRS"    /* RINEX marker name */
#define OSR_STAID   0                 /* RTCM 3 station ID */

#define MAXFILE     16                /* max number of input files */
#define OUT_FILE    "out.txt"         /* range corrections OSR default file */

#define SQR(x)      ((x)*(x))

#define UPDATE_OPT_ORBIT  (1 << 0)
#define UPDATE_OPT_CBIAS  (1 << 1)
#define UPDATE_OPT_PBIAS  (1 << 2)
#define UPDATE_OPT_L0BIAS (1 << 3)
#define UPDATE_OPT_CLCK   (1 << 4)
#define BUFFER_SIZE (16)
/* global variables ----------------------------------------------------------*/

static const char *usage[] = {
"usage: ssr2obs [options] file ...",
"",
"options: ([] as default)",
"  -ts y/m/d h:m:s   start time of OBS (GPST) []",
"  -te y/m/d h:m:s   end time of OBS   (GPST) [start time + 1h]",
"  -l6w week         specify GPS week corresponding to the start time of .l6 file [-ts time]",
"  -ti tint          time interval (s) [1]",
"  -k  file          configuration file []",
"  -o  file          output range corrections' OSR file [" OUT_FILE "]",
"  -dump             output each compact ssr subtype message",
"  -r                output RINEX3 OBS",
"  -b                output RTCM3 MSM",
"  -x                debug trace level (0: no trace) [0]",
"  file ...          QZSS L6 message file with extention '.l6' or '.L6' and",
"                    RINEX NAV files",
NULL
};

static rtcm_t rtcm_out = {0};         /* RTCM control struct */
static rnxopt_t rnx_opt = {0};        /* RINEX output options */

void output_osr_txt( FILE *fp, const obs_t *obs, const osrd_t *osr, const int n, double *pos);
extern int readotlgrid(const char *file, nav_t *nav);


/* set RINEX output options --------------------------------------------------*/
static void set_rnxopt(rnxopt_t *opt, char **infile, int n,
                       const prcopt_t *prcopt)
{
    static const int sys[] = {SYS_GPS, SYS_GLO, SYS_GAL, SYS_QZS};
    static const char *tobs[][OSR_NFREQ] = { /* supported obs types */
        {"1C", "2W", "2X", "5X"}, {"1C", "2P"}, {"1X", "5X"}, {"1C", "2X", "5X"}
    };
    gtime_t time = {0};
    int i, j, nobs;

    opt->rnxver = OSR_RNXVER;
    opt->navsys = prcopt->navsys;
    sprintf(opt->prog  , "%s %s", PROGNAME, PROG_VER);
    sprintf(opt->marker, "%s", OSR_MARKER);
    matcpy(opt->apppos, prcopt->ru, 3, 1);
    for (i = 0; i < n && i < MAXCOMMENT; i++) {
        sprintf(opt->comment[i], "infile: %-55.55s", infile[i]);
    }
    opt->tstart = opt->tend = time;
    
    for (i = 0; i < 4; i++) {
        if (!(sys[i] & opt->navsys)) continue;
        for (j = nobs = 0; j < OSR_NFREQ; j++) {
            if(tobs[i][j]==NULL) continue;
            sprintf(opt->tobs[i][nobs++], "C%s", tobs[i][j]);
            sprintf(opt->tobs[i][nobs++], "L%s", tobs[i][j]);
        }
        opt->nobs[i] = nobs;
    }
}

/* open output file ----------------------------------------------------------*/
static FILE *open_osr(const char *file, rtcm_t *rtcm, int mode, char **infile, int n,
                      const prcopt_t *prcopt)
{
    static const char osr_header[] =
        "msg,tow,sys,prn,pbias1,pbias2,pbias5,cbias1,cbias2,cbias5,trop,iono,antr1,antr2,antr5,relatv,wup1,wup2,wup5,compI1,compI2,compI5,compN,CPC1,CPC2,CPC5,PRC1,PRC2,PRC5,orb,clk,lat,lon,alt\n";
    FILE *fp=NULL;
    
    if (!(fp = fopen(file, mode == OSR_RTCM3 ? "wb" : "w"))) {
        fprintf(stderr, "Output file open error. %s\n", file);
        return NULL;
    }
    if (mode == OSR_RCORR) {
        fprintf(fp, "%s", osr_header);
    }
    else if (mode == OSR_RTCM3) {
        init_rtcm(&rtcm_out);
        rtcm_out.staid = OSR_STAID;
        matcpy(rtcm_out.sta.pos, prcopt->ru, 3, 1);
    }
    else if (mode == OSR_RINEX) {
        set_rnxopt(&rnx_opt, infile, n, prcopt);
        outrnxobsh(fp, &rnx_opt, &rtcm->nav);
    }
    return fp;
}

/* close output file ---------------------------------------------------------*/
static void close_osr(FILE *fp, rtcm_t *rtcm, int mode)
{

    if (mode == OSR_RTCM3) {
        free_rtcm(&rtcm_out);
    }
    else if (mode == OSR_RINEX) {
        rewind(fp);
        outrnxobsh(fp, &rnx_opt, &rtcm->nav);
    }
    fclose(fp);
}

/* write OSR to output file --------------------------------------------------*/
static void write_osr(FILE *fp, int mode, const obs_t *obs, const osrd_t *osr)
{
    static const int syss[] = {SYS_GPS, SYS_GLO, SYS_GAL, SYS_QZS};
    static const int types[] = {1074, 1084, 1094, 1114}; /* RTCM 3 MSM4 */
    double pos[3]={0.0, 0.0, 0.0};
    int i, j=0, sys;
    
    if (mode == OSR_RCORR) {
        output_osr_txt( fp, obs, osr, obs->n, pos);
    }
    else if (mode == OSR_RTCM3) {
        
        if (obs->n == 0) {
            return;
        }
        rtcm_out.time = obs->data[0].time;
        rtcm_out.obs.n = obs->n;
        
        for (i = sys = 0; i < obs->n; i++) {
            rtcm_out.obs.data[i] = obs->data[i];
            sys |= satsys(obs->data[i].sat, NULL);
        }
        for (i = 0; i < 4; i++) {
            if (sys & syss[i]) j = i;
        }
        for (i = 0; i < 4; i++) {
            if (!(sys & syss[i])) continue;
            
            if (gen_rtcm3(&rtcm_out, types[i], i < j)) {
                fwrite(rtcm_out.buff, rtcm_out.nbyte, 1, fp);
            }
        }
        if (gen_rtcm3(&rtcm_out, 1005, 0)) {
            fwrite(rtcm_out.buff, rtcm_out.nbyte, 1, fp);
        }
    }
    else if (mode == OSR_RINEX) {
        if (obs->n > 0) {
            outrnxobsb(fp, &rnx_opt, obs->data, obs->n, 0);
            
            if (!rnx_opt.tstart.time) {
                rnx_opt.tstart = obs->data[0].time;
            }
            rnx_opt.tend = obs->data[0].time;
        }
    }
}

/* calculate the actual distance between rcv and sat*/
static int actualdist(gtime_t time, obs_t *obs, nav_t *nav, const double *x)
{
    int i,j,n,sat,lsat[MAXSAT];
    double r,rr[3],dt,dt_p,*var,*dts,*rs,*e; 
    gtime_t tg;
    obsd_t *obsd = obs->data;
    int svh[MAXOBS];

    obs->n = 0;
    for (i=0;i<MAXOBS;i++) obsd[i].time = time;

    for (i=n=0;i<MAXSAT;i++) {
        if (!nav->ssr[i].t0[0].time||!nav->ssr[i].t0[1].time||
                nav->ssr[i].nsig==0) continue;
        lsat[n++]=i+1;
    }

    var=mat(1,n);dts=mat(2,n);rs=mat(6,n);  e = mat(n,3);
       
    for (i=0;i<3;i++) rr[i]=x[i];
    
    if (norm(rr,3)<=0.0) return -1;
    /* calculate code range */
    for (i=j=0;i<n;i++) {
        sat = lsat[i];
        dt = 0.08; dt_p = 0.0;
        while (1) {
            tg = timeadd(time,-dt);
            if(!ephpos(tg,time,sat,nav,nav->ssr[sat-1].iode,
                    rs+i*6,dts+2*i,var+i,svh)) {
                obsd[i].sat = 0;
                break;
            }
            if((r=geodist(rs+i*6,x,e+i*3))<=0.0) {
                obsd[i].sat = 0;
                break;
            }
            dt_p = dt;
            dt = r/CLIGHT;
            if (fabs(dt-dt_p)<1.0e-12) {
                obsd[i].time = time;
                obsd[i].sat = sat;
                obsd[i].P[0]=r+CLIGHT*(-dts[2*i]);
                break;
            }
        }
        j++;
    }
    free(var);  free(dts); free(e); free(rs);
    obs->n =n;
    return 0;
}


/* generate OSR --------------------------------------------------------------*/
static int gen_osr(gtime_t ts, gtime_t te, double ti, double latency, int mode,
                   prcopt_t *prcopt, filopt_t *fopt, solopt_t *solopt, 
                   char **infile, int n,
                   const char *outfile)
{
    static rtk_t rtk = {0};
    static rtcm_t rtcm = {0};
    static rtcm_t rtcm_ = {0};
    static obsd_t data[MAXOBS] = {0,};
    static osrd_t osr[MAXOBS] = {0,};
    FILE *fp_in, *fp_out, *dummy_fp[CSSR_TYPE_NUM];
    obs_t obs = {0};
    gtime_t time;
    char path[1024];
    int i,stat, nt = ti <= 0.0 ? 0 : (int)((timediff(te, ts)) / ti) + 1;
    int week_ref;
    int ret=0;

    /* read grid definition file */
    if (read_grid_def(fopt->grid)) {
        fprintf(stderr, "Grid file read error. %s\n", fopt->grid);
        showmsg("Grid file read error. %s\n", fopt->grid);
        return -1;
    }

    init_rtcm(&rtcm);
    init_rtcm(&rtcm_);

    /* read ocean tide loading parameters */
    if (prcopt->mode>PMODE_SINGLE&&fopt->blq) {
        if (prcopt->tidecorr == 3) {
            /* clasgrid */
            if(!readotlgrid(fopt->blq,&rtcm.nav)){
                freenav(&rtcm.nav,0xFF);
                return 0;
            }
        }
    }

    /* read RINEX NAV files */
    for (i = 0; i < n; i++) {
        stat = reppath(infile[i], path, ts, "", "");
        readrnx(path, 0, "", NULL, &rtcm.nav, NULL);
        
        if (stat) { /* next day with keywords */
            reppath(infile[i], path, timeadd(ts, 86400.0), "", "");
            readrnx(path, 0, "", NULL, &rtcm.nav, NULL);
        }
    }
    uniqnav(&rtcm.nav);
    
    if (rtcm.nav.n <= 0) {
        fprintf(stderr, "No navigation data.\n");
        return 0;
    }
    /* open QZSS L6 message file */
    if (!(fp_in = open_L6(infile, n))) {
        freenav(&rtcm.nav, 0xFF);
        return 0;
    }
    /* open output file */
    reppath(outfile, path, ts, "", "");
    
    if (!(fp_out = open_osr(path, &rtcm, mode, infile, n, prcopt))) {
        fclose(fp_in);
        freenav(&rtcm.nav, 0xFF);
        return 0;
    }
    
    rtkinit(&rtk, prcopt);
    obs.data = data;
    
    /* set reference week */
    if (prcopt->l6week != 0) {
        week_ref = prcopt->l6week;
    } else {
        time2gpst(ts,&week_ref);
    }
    for (i = 0; i < sizeof(rtcm_.week_ref) / sizeof(int); i++) {
        rtcm_.week_ref[i] = week_ref;
    }
    rtcm_.time=    ts;
    rtcm.nav.rtcmmode=rtcm_.nav.rtcmmode=RTCMMODE_CSSR;
    rtcm.mode=RTCMMODE_CSSR;

    for (i = 0; i < nt; i++) {
        time = timeadd(ts, ti * i); 
        rtk.sol.time= time;

        /* read compact SSR from L6 messsage file */
        while (rtcm_.subtype != 1 || timediff(rtcm_.time,time) <= -latency) {
            if ((ret=input_cssrf(&rtcm_,fp_in,dummy_fp))<-1) break;
        }

        if (rtcm_.nav.filreset == TRUE) {
                rtcm.nav.filreset = rtcm_.nav.filreset;
                rtcm_.nav.filreset = FALSE;
        }
        /* calculate the actual distance between rcv and sat*/
        if (actualdist(time,&obs,&rtcm.nav, prcopt->ru)<0) return -1;
        
        /* convert compact SSR to OSR */
        obs.n=ssr2osr(&rtk, obs.data, obs.n,&rtcm.nav,osr,mode);
        
        if (obs.n>0) {                
            /* write OSR */
            write_osr(fp_out, mode, &obs, osr);
        }
    }
    rtkfree(&rtk);
    free_rtcm(&rtcm);
    
    /* close output file */
    close_osr(fp_out, &rtcm, mode);
    
    fclose(fp_in);
    freenav(&rtcm.nav, 0xFF);
    
    return 0;
}

/* set processing options ----------------------------------------------------*/
static int set_prcopt(const char *file, prcopt_t *prcopt, solopt_t *solopt,
                            filopt_t *filopt, double *pos)
{
    prcopt->mode    = PMODE_SSR2OSR;
    prcopt->nf      = OSR_NFREQ;
    prcopt->navsys  = OSR_SYS;
    prcopt->elmin   = OSR_ELMASK * D2R;
    prcopt->tidecorr   = 1; /* tide correction */
    prcopt->posopt[2]  = 1; /* phase windup correction */
#if 0
    prcopt->cssropt[0] = 1; /* compensate time variation of ionosphere delay */
    prcopt->cssropt[1] = 1; /* shapiro time delay correction */
#endif
    solopt->timef=0;
    
    if (*file) {
        setsysopts(prcopt, solopt, filopt);
        if (!loadopts(file, sysopts)) {
            fprintf(stderr, "Configuration file read error. %s\n", file);
            return 0;
        }
        getsysopts(prcopt, solopt, filopt);
    }
    if (norm(pos, 3) > 0.0) {
        pos[0] *= D2R;
        pos[1] *= D2R;
        pos2ecef(pos, prcopt->ru);
    }
    if (norm(prcopt->ru, 3) <= 0.0) {
        fprintf(stderr, "No user position specified.\n");
        return 0;
    }
    return 1;
}


/* main ----------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts, te;
    double es[6] = {0,}, ee[6] = {0,}, ti = 1.0, latency = 0.0, pos[3] = {0,};
    char *conffile = "", *outfile = OUT_FILE;
    char *infile[MAXFILE];
    int i, j, n = 0, stat = 1, mode = OSR_RCORR;
    FILE *ofp[CSSR_TYPE_NUM];

    prcopt.mode  =PMODE_SSR2OSR;
    prcopt.navsys=SYS_GPS|SYS_QZS|SYS_GAL;
    prcopt.refpos=1;
    prcopt.glomodear=1;
    solopt.timef=0;
    sprintf(solopt.prog ,"%s ver.%s",PROGNAME,VER_RTKLIB);
    sprintf(filopt.trace,"%s.trace",PROGNAME);

    /* parse options */
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-ts") && i + 2 < argc) {
            sscanf(argv[++i], "%lf/%lf/%lf", es, es + 1,es + 2);
            sscanf(argv[++i], "%lf:%lf:%lf", es + 3, es + 4, es + 5);
            ts = epoch2time(es);
        }
        else if (!strcmp(argv[i], "-te") && i + 2 < argc) {
            sscanf(argv[++i], "%lf/%lf/%lf", ee, ee + 1, ee + 2);
            sscanf(argv[++i], "%lf:%lf:%lf", ee + 3, ee + 4, ee + 5);
            te=epoch2time(ee);
        }
        else if (!strcmp(argv[i],"-l6w")&&i+1<argc) prcopt.l6week=atof(argv[++i]);
        else if (!strcmp(argv[i], "-ti") && i + 1 < argc) ti = atof(argv[++i]);
        else if (!strcmp(argv[i], "-l") && i + 1 < argc) latency = atof(argv[++i]);
        else if (!strcmp(argv[i], "-p") && i + 1 < argc) sscanf(argv[++i], "%lf,%lf,%lf", pos, pos+1, pos+2);
        else if (!strcmp(argv[i], "-k") && i + 1 < argc) conffile = argv[++i];
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) outfile = argv[++i];
        else if (!strcmp(argv[i], "-x") && i + 1 < argc) solopt.trace = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-b")) mode = OSR_RTCM3;
        else if (!strcmp(argv[i], "-r")) mode = OSR_RINEX;
        else if (!strcmp(argv[i],"-dump")) mode=OSR_L6DUMP;
        else if (!strncmp(argv[i], "-", 1)) {
            for (j = 0; usage[j]; j++) fprintf(stderr, "%s\n", usage[j]);
            return 0;
        }
        else if (n < MAXFILE) infile[n++] = argv[i];
    }
    if (norm(es, 3) <= 0.0 && mode != OSR_L6DUMP) {
        fprintf(stderr, "No start time specified.\n");
        return -1;
    }
    if (norm(ee, 3) <= 0.0) {
        te = timeadd(epoch2time(es), 3599.0);
    }
    if (!set_prcopt(conffile, &prcopt,&solopt,&filopt, pos)) return -1;
    if (solopt.trace > 0) {
        traceopen(filopt.trace);
        tracelevel(solopt.trace);
    }
    
    if (mode == OSR_L6DUMP) { 
        /* parse cssr */
        if (open_outputfiles(ofp) == -1) {
            fprintf(stderr, "Can't open output files.\n");
            return -1;
        }
        stat = dumpcssr(infile, n, ofp, &filopt);
        close_outputfiles(ofp);
    } else {
        /* generate OSR */
        stat = gen_osr(ts, te, ti, latency, mode, &prcopt, &filopt, 
                    &solopt, infile, n, outfile);
    }
    traceclose();

    return stat;
}
