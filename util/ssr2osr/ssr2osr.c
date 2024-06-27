/*------------------------------------------------------------------------------
* ssr2osr.c : read rinex obs/nav files, QZSS L6 Message(.l6) 
*             and output individual corrections of observation space representation(OSR).
*             Derived from rnx2rktp.
*
*          Copyright (C) 2007- by T.TAKASU, All rights reserved.
*          Copyright (C) 2015- by Mitsubishi Electric Corp., All rights reserved.
*
* version : $Revision:  1.0
* history   2018/3/29   1.0 new
*
* Reference(rtklib)
* history : 2007/01/16  1.0 new
*           2007/03/15  1.1 add library mode
*           2007/05/08  1.2 separate from postpos.c
*           2009/01/20  1.3 support rtklib 2.2.0 api
*           2009/12/12  1.4 support glonass
*                           add option -h, -a, -l, -x
*           2010/01/28  1.5 add option -k
*           2010/08/12  1.6 add option -y implementation (2.4.0_p1)
*           2014/01/27  1.7 fix bug on default output time format
*           2015/05/15  1.8 -r or -l options for fixed or ppp-fixed mode
*           2015/06/12  1.9 output patch level in header
*           2016/09/07  1.10 add option -sys
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"
#include "cssr.h"

#define PROGNAME    "ssr2osr"           /* program name */
#define MAXFILE     8                   /* max number of input files */

/* help text -----------------------------------------------------------------*/
static const char *help[]={
"",
" usage: ssr2osr [option]... file file [...]",
"",
" Read RINEX OBS/NAV, QZSS L6 MESSAGE(.l6) files and output",
" individual corrections of observation space representation(OSR).",
" -?              print help",
" -k   file       input options from configuration file [off]",
" -ti  tint       time interval (sec) [all]",
" -ts  ds ts      start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te  de te      end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -l6w week       specify GPS week corresponding to the start time of .l6 file [obs start time]",
" -o   file       set output file [NMEA-GGA]",
" -dump           output each compact ssr subtype message",
" -x   level      debug trace level (0:off) [0]"
};

/* show message --------------------------------------------------------------*/
extern int showmsg(char *format, ...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}
/* ssr2osr main -------------------------------------------------------------*/
int main(int argc, char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts={0},te={0};
    double tint=1.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59},pos[3];
    int i,j,n,ret=0;
    char *infile[MAXFILE],*outfile="";
    FILE *ofp[CSSR_TYPE_NUM];
    
    prcopt.mode  =PMODE_SSR2OSR;
    prcopt.navsys=SYS_GPS|SYS_QZS;
    prcopt.refpos=1;
    prcopt.glomodear=1;
    solopt.timef=0;
    sprintf(solopt.prog ,"%s ver.%s",PROGNAME,VER_RTKLIB);
    sprintf(filopt.trace,"%s.trace",PROGNAME);
    
    /* load options from configuration file */
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-k")&&i+1<argc) {
            resetsysopts();
            if (!loadopts(argv[++i],sysopts)) return -1;
            getsysopts(&prcopt,&solopt,&filopt);
        }
    }
    for (i=1,n=0;i<argc;i++) {
        if      (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-ts")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",es,es+1,es+2);
            sscanf(argv[++i],"%lf:%lf:%lf",es+3,es+4,es+5);
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",ee,ee+1,ee+2);
            sscanf(argv[++i],"%lf:%lf:%lf",ee+3,ee+4,ee+5);
            te=epoch2time(ee);
        }
        else if (!strcmp(argv[i],"-l6w")&&i+1<argc) prcopt.l6week=atof(argv[++i]);
        else if (!strcmp(argv[i],"-ti")&&i+1<argc) tint=atof(argv[++i]);
        else if (!strcmp(argv[i],"-k")&&i+1<argc) {++i; continue;}
        else if (!strcmp(argv[i],"-x")&&i+1<argc) solopt.trace=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-dump")) {
            prcopt.l6dump=1;
        }
        else if (*argv[i]=='-') printhelp();
        else if (n<MAXFILE) infile[n++]=argv[i];
    }
    if (n<=0) {
        showmsg("error : no input file");
        return -2;
    }

    if (prcopt.l6dump) {
        if (open_outputfiles(ofp) == -1) {
            fprintf(stderr, "Can't open output file.\n");
            return -1;
        }
        /* parse cssr */
        ret=dumpcssr(infile, n, ofp, &filopt);
        close_outputfiles(ofp);
    } else {
        /* generate OSR */
        ret=postpos(ts,te,tint,0.0,&prcopt,&solopt,&filopt,infile,n,outfile,"","");
        if (!ret) fprintf(stderr,"%40s\r","");
    }
    
    return ret;
}
