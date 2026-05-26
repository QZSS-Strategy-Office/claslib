/*------------------------------------------------------------------------------
* ppprtk.c : precise point positioning real time kinematic
*
*          Copyright (C) 2007- by T.TAKASU, All rights reserved.
*          Copyright (C) 2015- by Mitsubishi Electric Corp., All rights reserved.
*
* options : -DIERS_MODEL use IERS tide model
*
* references :
*     [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*         2003, November 2003
*     [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*         celestial reference frames, Geophys. Res. Let., 1997
*     [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*
* version : $Revision:$ $Date:$
* history   2018/3/29   1.0 new as claslib
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "cssr2osr.h"

#define INTERPMODE 1
#define MAXREF 2                    /* max interation of searching reference sat */
#define SQR(x)      ((x)*(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))

#define AS2R        (D2R/3600.0)    /* arc sec to radian */
#define GME         3.986004415E+14 /* earth gravitational constant */
#define GMS         1.327124E+20    /* sun gravitational constant */
#define GMM         4.902801E+12    /* moon gravitational constant */
                                    
                                    /* initial variances */
#define VAR_POS     SQR(30.0)       /* initial variance of receiver pos ((m)^2) */
#define VAR_VEL     SQR(1.0)        /* initial variance of receiver acc ((m/s)^2) */
#define VAR_ACC     SQR(1.0)        /* initial variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(100.0)      /* receiver clock (m^2) */
#define VAR_ZTD     SQR(  0.3)      /*   ztd (m^2) */
#define VAR_GRA     SQR(0.001)      /*   gradient (m^2) */
#define INIT_ZWD    0.001           /* initial zwd (m) */

#define GAP_RESION  120      /* gap to reset ionosphere parameters (epochs) */
#define BASELEN  1000.0      /* baseline length(m):1km=>5mm,2km->10mm,10km->50mm */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */

#define MAX_SAT 14                  /* number of satellites in SSR */
#define MAX_NGRID 4                 /* number of gridded points for interpolation */


#define MAXPBCORSSR 20.0            /* max phase bias correction of ssr (m) */
#define CSSRINVALID -10000          /* invalid value */

static const double sfi[4]={FREQ2/FREQ1,FREQ1/FREQ2,1.0-SQR(FREQ1/FREQ2),FREQ2/FREQ1-FREQ1/FREQ2};

/*  The F inverse distribution function */
static const double qf[6][60] = {
 /* significance level, alpha=0.001(0.1%) */
    {405284.07,999.00,141.11,53.44,29.75,20.03,15.02,12.05,10.11,8.75,
     7.76,7.00,6.41,5.93,5.54,5.20,4.92,4.68,4.47,4.29,
     4.13,3.98,3.85,3.74,3.63,3.53,3.44,3.36,3.29,3.22,
     3.15,3.09,3.04,2.98,2.93,2.89,2.84,2.80,2.76,2.73,
     2.69,2.66,2.63,2.60,2.57,2.54,2.51,2.49,2.46,2.44,
     2.42,2.40,2.38,2.36,2.34,2.32,2.30,2.28,2.27,2.25},
 /* significance level, alpha=0.005(0.5%) */
    {16210.72,199.00,47.47,23.15,14.94,11.07,8.89,7.50,6.54,5.85,
    5.32,4.91,4.57,4.30,4.07,3.87,3.71,3.56,3.43,3.32,
    3.22,3.12,3.04,2.97,2.90,2.84,2.78,2.72,2.67,2.63,
    2.58,2.54,2.51,2.47,2.44,2.41,2.38,2.35,2.32,2.30,
    2.27,2.25,2.23,2.21,2.19,2.17,2.15,2.13,2.11,2.10,
    2.08,2.07,2.05,2.04,2.02,2.01,2.00,1.99,1.97,1.96},
 /* significance level, alpha=0.010(1.0%) */
    {4052.18,99.00,29.46,15.98,10.97,8.47,6.99,6.03,5.35,4.85,
    4.46,4.16,3.91,3.70,3.52,3.37,3.24,3.13,3.03,2.94,
    2.86,2.78,2.72,2.66,2.60,2.55,2.51,2.46,2.42,2.39,
    2.35,2.32,2.29,2.26,2.23,2.21,2.18,2.16,2.14,2.11,
    2.09,2.08,2.06,2.04,2.02,2.01,1.99,1.98,1.96,1.95,
    1.94,1.92,1.91,1.90,1.89,1.88,1.87,1.86,1.85,1.84},
 /* significance level, alpha=0.050(5.0%) */
    {161.45,19.00,9.28,6.39,5.05,4.28,3.79,3.44,3.18,2.98,
     2.82,2.69,2.58,2.48,2.40,2.33,2.27,2.22,2.17,2.12,
     2.08,2.05,2.01,1.98,1.96,1.93,1.90,1.88,1.86,1.84,
     1.82,1.80,1.79,1.77,1.76,1.74,1.73,1.72,1.70,1.69,
     1.68,1.67,1.66,1.65,1.64,1.63,1.62,1.62,1.61,1.60,
     1.59,1.58,1.58,1.57,1.56,1.56,1.55,1.55,1.54,1.53},
 /* significance level, alpha=0.100(10.0%) */       
    {39.86,9.00,5.39,4.11,3.45,3.05,2.78,2.59,2.44,2.32,
    2.23,2.15,2.08,2.02,1.97,1.93,1.89,1.85,1.82,1.79,
    1.77,1.74,1.72,1.70,1.68,1.67,1.65,1.63,1.62,1.61,
    1.59,1.58,1.57,1.56,1.55,1.54,1.53,1.52,1.51,1.51,
    1.50,1.49,1.48,1.48,1.47,1.46,1.46,1.45,1.45,1.44,
    1.44,1.43,1.43,1.42,1.42,1.41,1.41,1.40,1.40,1.40},
 /* significance level, alpha=0.200(20.0%) */   
    {9.47,4.00,2.94,2.48,2.23,2.06,1.94,1.86,1.79,1.73,
    1.69,1.65,1.61,1.58,1.56,1.54,1.52,1.50,1.48,1.47,
    1.45,1.44,1.43,1.42,1.41,1.40,1.39,1.38,1.37,1.36,
    1.36,1.35,1.34,1.34,1.33,1.33,1.32,1.32,1.31,1.31,
    1.30,1.30,1.30,1.29,1.29,1.28,1.28,1.28,1.27,1.27,
    1.27,1.26,1.26,1.26,1.26,1.25,1.25,1.25,1.25,1.24}
};

extern FILE *fp_osr;
#ifdef CSSR2OSR_VRS
static unsigned char obsfreqs[]={ /* 1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,7:L3 */
    
    0, 1, 1, 1, 1,  1, 1, 1, 1, 1, /*  0- 9 */
    1, 1, 1, 1, 2,  2, 2, 2, 2, 2, /* 10-19 */
    2, 2, 2, 2, 3,  3, 3, 5, 5, 5, /* 20-29 */
    4, 4, 4, 4, 4,  4, 4, 6, 6, 6, /* 30-39 */
    2, 2, 4, 4, 3,  3, 3, 1, 1, 0  /* 40-49 */
};
#endif


/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int type, const prcopt_t *opt)
{
    double a,b,a2,b2,c,d,fact=1.0;
    double sinel=sin(el);
    int i=sys==SYS_GLO?1:(sys==SYS_GAL?2:0);

    c=opt->ionoopt==IONOOPT_EST?opt->err[5]:0.0;
    d=opt->tropopt>=TROPOPT_EST?opt->err[6]:0.0;

    /* extended error model */
    if (type==1&&opt->exterr.ena[0]) { /* code */
        a=opt->exterr.cerr[i][0];
        b=opt->exterr.cerr[i][1];
        if (opt->ionoopt==IONOOPT_IFLC) {
            a2=opt->exterr.cerr[i][2];
            b2=opt->exterr.cerr[i][3];
            a=sqrt(SQR(2.55)*a*a+SQR(1.55)*a2*a2);
            b=sqrt(SQR(2.55)*b*b+SQR(1.55)*b2*b2);
        }
    }
    else if (type==0&&opt->exterr.ena[1]) { /* phase */
        a=opt->exterr.perr[i][0];
        b=opt->exterr.perr[i][1];
        if (opt->ionoopt==IONOOPT_IFLC) {
            a2=opt->exterr.perr[i][2];
            b2=opt->exterr.perr[i][3];
            a=sqrt(SQR(2.55)*a*a+SQR(1.55)*a2*a2);
            b=sqrt(SQR(2.55)*b*b+SQR(1.55)*b2*b2);
        }
        /*c=d=0.0;*/

    }
    else if (type==2) { /* iono */
        a=b=d=0.0;
    }
    else if (type==3) { /* trop */
        a=b=c=0.0;
    }
    else { /* normal error model */
        if (type==1) fact*=opt->eratio[0];
        fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
        if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;
        a=fact*opt->err[1];
        b=fact*opt->err[2];
        c=d=0.0;
    }
    return a*a+b*b/sinel/sinel+c*c+d*d;
}

/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}
extern int selfreqpair(const int sat, const prcopt_t *opt,const obsd_t *obs)
{
    int optf=opt->posopt[10];
    if (NFREQ==1||optf==POSL1) {
         return 0;
    }else {
        if ((satsys(sat,NULL))&(SYS_GAL)) {/* GAL */
            if(obs->L[2]!=0.0||obs->P[2]!=0.0) return 2; 
            return 0; 
        }
        if (optf==POSL1L2L5) return 1+2;
        if (optf==POSL1L5) return 2;
        if (optf==POSL1L5_L2&&obs->L[2]!=0.0&&obs->P[2]!=0.0) return 2;
        if (optf==POSL1L2) {
             if (obs->L[1]!=0.0||obs->P[1]!=0.0) return 1;
             return 0;
        }
        return 1;
    }
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t *rtk, const obsd_t *obs, int n,  double tt)
{
    double *F,*FP,*xp,pos[3],Q[9]={0},Qv[9],var=0.0;
    int i, j, na,nb,nx,flag_init=0;
    double *Fx,*Pxx,*Pxi,*Pxb,*Pxi_,*Pxb_,*xp_,ki;
    int ni;
    
    trace(3,"udpos_ppp:\n");

    nx=rtk->nx;

    /* fixed mode */
    if (rtk->opt.mode==PMODE_PPP_FIXED) {
        for (i=0;i<3;i++) initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }

    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
        flag_init=1;
    }

    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*nx];
    var/=3.0;

    if (var>VAR_POS) {
        /* reset position with large variance */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }

    /* state transition of position/velocity/acceleration */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        F=eye(nx);
        
        /* dynamics on */
        if (rtk->opt.dynamics) {
            for (i=0;i<6;i++) F[i+(i+3)*nx]=tt;
            for (i=0;i<3;i++) F[i+(i+6)*nx]=SQR(tt)/2.0;
        }

        na=NP(&rtk->opt);ni=MAXSAT;nb=nx-na-ni;
        ki=rtk->opt.beta>0.0?exp(-tt/rtk->opt.beta):1.0;
        Fx=mat(na,na);xp=mat(na,1);xp_=mat(na,1);FP=mat(na,na);
        Pxx=mat(na,na);Pxi=mat(na,ni);Pxi_=mat(na,ni);Pxb=mat(na,nb);Pxb_=mat(na,nb); 

        for (i=0;i<na;i++) for (j=0;j<na;j++) Pxx[i+j*na]=rtk->P[i+j*nx];
        for (i=0;i<na;i++) for (j=0;j<ni;j++) Pxi_[i+j*na]=0.5*(rtk->P[i+(na+j)*nx]+rtk->P[na+j+(i)*nx]);
        for (i=0;i<na;i++) for (j=0;j<nb;j++) Pxb_[i+j*na]=0.5*(rtk->P[i+(na+ni+j)*nx]+rtk->P[na+ni+j+(i)*nx]);
        for (i=0;i<na;i++) for (j=0;j<na;j++) Fx[i+j*na]=F[i+j*nx];
        for (i=0;i<na;i++) xp_[i]=rtk->x[i];

        /* x=F*x */
        matmul("NN",na,1,na,1.0,Fx,xp_,0.0,xp);
        for (i=0;i<na;i++) rtk->x[i]=xp[i];
        for (i=0;i<ni;i++) rtk->x[i+na]*=ki;

        /* Pxx=Fx*Pxx*Fx */
        matmul("NN",na,na,na,1.0,Fx,Pxx,0.0,FP);
        matmul("NT",na,na,na,1.0,FP,Fx,0.0,Pxx);
        for (i=0;i<na;i++) for (j=0;j<na;j++) rtk->P[i+j*nx]=Pxx[i+j*na];
        /* Pxi=k*Fx*Pxi*/
        matmul("NN",na,ni,na,ki,Fx,Pxi_,0.0,Pxi);
        for (i=0;i<na;i++) for (j=0;j<ni;j++) {
            rtk->P[i+(na+j)*nx]=rtk->P[(na+j)+i*nx]=Pxi[i+j*na];
        }
        /* Pxb=Fx*Pxb*/
        matmul("NN",na,nb,na,1.0,Fx,Pxb_,0.0,Pxb);
        for (i=0;i<na;i++) for (j=0;j<nb;j++) {
            rtk->P[i+(na+ni+j)*nx]=rtk->P[(na+ni+j)+i*nx]=Pxb[i+j*na];
        }

        for (i=0;i<ni;i++) for (j=0;j<ni;j++) rtk->P[(i+na)+(j+na)*nx]*=(ki*ki);
        for (i=0;i<ni;i++) for (j=0;j<nb;j++) rtk->P[(i+na)+(j+na+ni)*nx]*=ki;
        for (i=0;i<nb;i++) for (j=0;j<ni;j++) rtk->P[(i+na+ni)+(j+na)*nx]*=ki;

        free(Pxx);free(Pxi);free(Pxb);free(Pxi_);free(Pxb_);free(Fx);free(xp_);
        free(F); free(FP); free(xp);
    } else if (rtk->opt.dynamics) {
        /* calculate state transition of position/velocity/acceleration */
        /* x=Fx, P=FPF' */
        F=eye(9);  xp=mat(9,1); FP=mat(nx,nx);
        /* generate F matrix to update position and velocity */
        for (i=0;i<6;i++) F[i+(i+3)*9]=tt;
        for (i=0;i<3;i++) F[i+(i+6)*9]=SQR(tt)/2.0;

        /* x=F*x, only calculate pos/vel/acc states to save time, the rest are unchanged */
        matmul("NN",9,1,9,1.0,F,rtk->x,0.0,xp);
        matcpy(rtk->x,xp,9,1);
        /* P=F*P */
        matcpy(FP,rtk->P,nx,nx);
        for (j=0;j<nx;j++) for (i=0;i<6;i++) FP[i+j*nx]+=rtk->P[i+3+j*nx]*tt;
        /* P=FP*F'*/
        matcpy(rtk->P,FP,nx,nx);
        for (j=0;j<nx;j++) for (i=0;i<6;i++) rtk->P[j+i*nx]+=rtk->P[j+(i+3)*nx]*tt;
        free(F); free(FP); free(xp);
    }

    if (rtk->opt.dynamics) {
        /* process noise added to only acceleration */
        Q[0]=Q[4]=SQR(rtk->opt.prn[3])*fabs(tt); Q[8]=SQR(rtk->opt.prn[4])*fabs(tt);
        ecef2pos(rtk->x,pos);
        covecef(pos,Q,Qv);
        for (i=0;i<3;i++) for (j=0;j<3;j++) {
            rtk->P[i+6+(j+6)*nx]+=Qv[i+j*3];
        }
        if (rtk->opt.prnadpt&&flag_init==0) {
            for (i=0;i<9;i++) for (j=0;j<9;j++) {
                rtk->P[i+j*nx]+=rtk->Q[i+j*nx]*fabs(tt);
            }
        }
    } else {
        /* process noise added to only position */
        Q[0]=Q[4]=SQR(rtk->opt.prn[5])*fabs(tt); Q[8]=SQR(rtk->opt.prn[6])*fabs(tt);
        ecef2pos(rtk->x,pos);
        covecef(pos,Q,Qv);
        if (rtk->opt.prnadpt==0||flag_init) {
            for (i=0;i<3;i++) for (j=0;j<3;j++) {
                rtk->P[i+j*nx]+=Qv[i+j*3];
            }
        } else {
            for (i=0;i<3;i++) for (j=0;j<3;j++) {
                rtk->P[i+j*nx]+=rtk->Q[i+j*nx]*fabs(tt);
            }
        }
    }
}

/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop(rtk_t *rtk, double tt)
{
    int i,j,k;

    trace(3,"udtrop  : tt=%.1f\n",tt);

    for (i=0;i<1;i++) {
        j=ITT(&rtk->opt);

        if (rtk->x[j]==0.0) {
            initx(rtk,INIT_ZWD,SQR(rtk->opt.std[2]),j); /* initial zwd */
            rtk->Q[j+j*rtk->nx]=SQR(rtk->opt.prn[2])*tt;
            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) initx(rtk,1E-6,VAR_GRA,++j);
            }
        }
        else {
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2])*tt;
            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) {
                    rtk->P[++j*(1+rtk->nx)]+=SQR(rtk->opt.prn[2]*0.3)*fabs(rtk->tt);
                }
            }
        }
    }
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udion(rtk_t *rtk, double tt, double bl, const obsd_t *obs, int ns)
{
    double el,fact,qi;
    int i,j,k,f,m;

    trace(3,"udion  : tt=%.1f bl=%.0f ns=%d\n",tt,bl,ns);

    for (i=1;i<=MAXSAT;i++) {
        j=II(i,&rtk->opt);
        if (rtk->x[j]==0.0) continue;
        for (k=0;k<ns;k++) {
            if (obs[k].sat==i) break;
        }
        m=selfreqpair(i,&rtk->opt,obs+k);
        for (f=0;f<rtk->opt.nf;f++) {
            if (f>0&&!(f&m)) continue;
            if ((rtk->ssat[i-1].outc[f]>(unsigned int)rtk->opt.maxout)) {
                trace(3,"reset iono estimation, sat=%2d outc=%6d \n",i,
                      rtk->ssat[i-1].outc[f]);
                rtk->x[j]=0.0;
            }
        }
    }
    for (i=0;i<ns;i++) {
        j=II(obs[i].sat,&rtk->opt);
        if (rtk->x[j]==0.0) {
            initx(rtk,1E-6,SQR(rtk->opt.std[1]),j);
            rtk->Q[j+j*rtk->nx]=0.0;
        }
        else {
            /* elevation dependent factor of process noise */
            el=rtk->ssat[obs[i].sat-1].azel[1];
            fact=cos(el);
            if (rtk->opt.ionoopt==IONOOPT_EST_ADPT) {
                if (rtk->Q[j+j*rtk->nx]==0.0) {
                    rtk->Q[j+j*rtk->nx]=SQR(rtk->opt.prn[1]*fact);
                } else {
                    if (rtk->Q[j+j*rtk->nx]>SQR(rtk->opt.prn[7])) {
                        rtk->Q[j+j*rtk->nx]=SQR(rtk->opt.prn[7]);
                    } else if (rtk->Q[j+j*rtk->nx]<SQR(rtk->opt.prn[1])) {
                        rtk->Q[j+j*rtk->nx]=SQR(rtk->opt.prn[1]);
                    }
                }
                qi= rtk->Q[j+j*rtk->nx];
            } else {
                qi=SQR(rtk->opt.prn[1]*fact);
                rtk->Q[j+j*rtk->nx]=qi;
            }
            
            rtk->P[j+j*rtk->nx]+=qi*tt;
        }
    }
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int n)
{
    int i,j;

    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<rtk->opt.nf;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;

        trace(2,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);

        rtk->ssat[obs[i].sat-1].slip[j]=1;
    }
}

/* L1/L2 geometry-free phase measurement -------------------------------------*/
static double gfmeasL1L2(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    if (lam[0]==0.0||lam[1]==0.0||obs->L[0]==0.0||obs->L[1]==0.0) return 0.0;
    return lam[0]*obs->L[0]-lam[1]*obs->L[1];
}
/* L1/L5 geometry-free phase measurement -------------------------------------*/
static double gfmeasL1L5(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    if (lam[0]==0.0||lam[2]==0.0||obs->L[0]==0.0||obs->L[2]==0.0) return 0.0;
    return lam[0]*obs->L[0]-lam[2]*obs->L[2];
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double g0,g1;
    int i;
    
    trace(3,"detslip_gf_L1_L2\n");
    
    /* slip detection for L1-L2 gf  */
    for (i=0;i<n&&i<MAXOBS;i++) {
        if ((g1=gfmeasL1L2(obs+i,nav))==0.0) continue;
        g0=rtk->ssat[obs[i].sat-1].gf;
        rtk->ssat[obs[i].sat-1].gf=g1;
        if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
            trace(2,"detslip_gf_L1_L2: slip detected sat=%2d gf=%8.3f->%8.3f\n",
                  obs[i].sat,g0,g1);
            rtk->ssat[obs[i].sat-1].slip[0]|=1;
            rtk->ssat[obs[i].sat-1].slip[1]|=1;
        }
    }
    
    trace(3,"detslip_gf_L1_L5\n");
    
    /* slip detection for L1-L5 gf  */
    if (rtk->opt.nf<=2) return;
    for (i=0;i<n&&i<MAXOBS;i++) {
        if ((g1=gfmeasL1L5(obs+i,nav))==0.0) continue;
        g0=rtk->ssat[obs[i].sat-1].gf2;
        rtk->ssat[obs[i].sat-1].gf2=g1;
        if (g0!=0.0&&fabs(g1-g0)>rtk->opt.thresslip) {
            trace(2,"detslip_gf_L1_L5: slip detected sat=%2d gf=%8.3f->%8.3f\n",
                  obs[i].sat,g0,g1);
            rtk->ssat[obs[i].sat-1].slip[0]|=1;
            rtk->ssat[obs[i].sat-1].slip[2]|=1;
        }
    }
}

/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t *rtk, double tt, const obsd_t *obs, int n,
                       const nav_t *nav)
{
    double offset,lami,cp,pr;
    double *bias;
    int i,j,f,slip,nf=NF(&rtk->opt),reset,isat;

    trace(3,"udbias  : tt=%.1f n=%d\n",tt,n);

    for (i=0;i<MAXSAT;i++) for (j=0;j<rtk->opt.nf;j++) {
        rtk->ssat[i].slip[j]=0;
        rtk->ssat[i].pbreset[j]=0;
    }

    /* detect cycle slip by LLI */
    detslp_ll(rtk,obs,n);

    /* detect cycle slip by geometry-free phase jump */
    detslp_gf(rtk,obs,n,nav);

    for (f=0;f<nf;f++) {
        /* reset phase-bias if expire obs outage counter */
        for (i=0;i<MAXSAT;i++) {

            reset=++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout;

            if (rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i+1,f,&rtk->opt)]!=0.0) {
                initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));
            }
            else if (reset&&rtk->x[IB(i+1,f,&rtk->opt)]!=0.0) {
                /* reset udbias */
                initx(rtk,0.0,0.0,IB(i+1,f,&rtk->opt));
                trace(3,"udbias : obs outage counter overflow (sat=%3d L%d n=%5d)\n",
                      i+1,f+1,rtk->ssat[i].outc[f]);
            }
            if (rtk->opt.modear!=ARMODE_INST&&reset) {
                rtk->ssat[i].lock[f]=-rtk->opt.minlock;
            }
        }

        /* reset phase-bias if detecting cycle slip */
        for (i=0;i<n&&i<MAXOBS;i++) {
            j = IB(obs[i].sat,f,&rtk->opt);
            rtk->P[j+j*rtk->nx] += rtk->opt.prn[0]*rtk->opt.prn[0]*tt;
            slip = rtk->ssat[obs[i].sat-1].slip[f];
            if (rtk->opt.ionoopt==IONOOPT_IFLC) slip|=rtk->ssat[obs[i].sat-1].slip[1];
            if (rtk->opt.modear==ARMODE_INST||!(slip&1)) continue;
            if (!(slip&1)) continue;
            rtk->x[j] = 0.0;
            rtk->ssat[obs[i].sat-1].lock[f] = -rtk->opt.minlock;
        }
        bias=zeros(n,1);

        /* estimate approximate phase-bias by phase - code */
        for (i=j=0,offset=0.0;i<n;i++) {
            cp = obs[i].L[f];
            pr = obs[i].P[f];
            lami = nav->lam[obs[i].sat-1][f];
            if (cp==0.0||pr==0.0||lami==0.0) continue;

            bias[i]=cp*lami-pr;

            /* offset = sum of (bias - phase-bias) for all valid sats in meters */
            if (rtk->x[IB(obs[i].sat,f,&rtk->opt)]!=0.0) {
                lami=nav->lam[obs[i].sat-1][f];
                offset+=bias[i]-rtk->x[IB(obs[i].sat,f,&rtk->opt)]*lami;
                j++;
            }
        }

        rtk->com_bias=j>0?offset/j:0; /* save offset for initialization below */
        /* set initial states of phase-bias */
        for (i=0;i<n;i++) {
            isat=obs[i].sat;
            if (bias[i]==0.0||rtk->x[IB(isat,f,&rtk->opt)]!=0.0) continue;
            trace(3,"set initial states of phase-bias, sat=%2d f=%2d\n",obs[i].sat,f+1);
            lami=nav->lam[isat-1][f];
            initx(rtk,(bias[i]-rtk->com_bias)/lami,SQR(rtk->opt.std[0]),IB(obs[i].sat,f,&rtk->opt));
            rtk->ssat[isat-1].lock[f]=-rtk->opt.minlock;
            rtk->ssat[isat-1].pbreset[f]=1;
        }

        free(bias);
    }
}
/* temporal update of states --------------------------------------------------*/
extern void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    double tt=fabs(rtk->tt),bl;

    trace(3,"udstate_ppp: n=%d tt=%f \n",n,tt);

    /* temporal update of position */
    udpos_ppp(rtk,obs,n,tt);

    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt>=TROPOPT_EST) {
        udtrop(rtk,tt);
    }

    /* temporal update of phase-bias */
    udbias_ppp(rtk,tt,obs,n,nav);

    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        bl=BASELEN; /* baseline length*/
        udion(rtk,tt,bl,obs,n);
    }
}


/* precise tropspheric model -------------------------------------------------*/
static double prectrop2(gtime_t time, const double *pos, int r,
                        const double *azel, const prcopt_t *opt,
                        const double *x, double *dtdx)
{
    double m_w=0.0,cotz,grad_n,grad_e;
    int i=ITT(opt);
    
    /* wet mapping function */
    tropmapf(time,pos,azel,&m_w);
    
    if (opt->tropopt>=TROPOPT_ESTG&&azel[1]>0.0) {
        
        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[i+1]+grad_e*x[i+2];
        dtdx[1]=grad_n*x[i];
        dtdx[2]=grad_e*x[i];
    }
    else dtdx[1]=dtdx[2]=0.0;
    dtdx[0]=m_w;
    return m_w*x[i];
}

static int validobs(int i,int f, int nf, double *y)
{
    return y[f+i*nf*2] != 0.0 && (f<nf || (y[f-nf+i*nf*2] != 0.0));
}

static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                  int nv, double *R)
{
    int i,j,k=0,b;

    for (i=0;i<nv*nv;i++) R[i] = 0.0;
    for (b=0;b<n;k+=nb[b++]) {
        for (i=0;i<nb[b];i++) {
            for (j=0;j<nb[b];j++) {
                R[k+i+(k+j)*nv] = Ri[k+i]+((i==j) ? Rj[k+i] : 0.0);
            }
        }
    }
}
/* test navi system (m=0:gps/qzs/sbs,1:glo,2:gal,3:bds,4:qzs,5:gps(l2c)) ------*/
static int test_sys(int sys, int m, int qzsmodear, int code)
{
    switch (sys) {
        case SYS_GPS: return m==(code == CODE_L2X ? 5 : 0);
        case SYS_QZS: return m==(qzsmodear == 2 ? 0 : 4);
        case SYS_SBS: return m==0;
        case SYS_GLO: return m==1;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
    }
    return 0;
}
/* single-differenced phase/code residuals -----------------------------------*/
extern int ddres(rtk_t *rtk, const nav_t *nav, const double *x,
                 const double *P, const obsd_t *obs, double *y, double *e,
                 double *azel, int n, double *v, double *H, double *R,
                 int *vflg, int niter)
{
    prcopt_t *opt=&rtk->opt;
    double pos[3],lami,lamj,*Ri,*Rj,*Hi=NULL;
    double *tropu,*im,*dtdxu,didxi=0.0,didxj=0.0,fi,fj;
    static int refsat[NFREQ*2*MAXREF*6+1]={0};
    int i,j,k,m,l,f,ff,nv=0,nb[NFREQ*6*2+1]={0},b=0,sati,satj,sysi,sysj,nf=NF(opt);
    int flg=0,code;

    trace(3,"ddres   :niter=%2d nx=%d n=%d\n",niter,rtk->nx,n);

    im=mat(n,1);tropu=mat(n,1);dtdxu=mat(n,3);Ri=mat(n*nf*2*2,1);Rj=mat(n*nf*2*2,1);
    ecef2pos(x,pos);

    for (i=0;i<MAXSAT;i++) {
        for (j=0;j<NFREQ;j++) rtk->ssat[i].resp[j]=rtk->ssat[i].resc[j]= 0.0;
    }
    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<n;i++) {
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=ionmapf(pos,azel+i*2);
        }
        if (opt->tropopt>=TROPOPT_EST) {
            tropu[i]=prectrop2(rtk->sol.time,pos,0,azel+i*2,opt,x,dtdxu+i*3);
        }
    }

    for (m=0;m<6;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds,4:qzs */
        for (j=l=0;j<n;j++) {
            sati=obs[j].sat;
            code=rtk->ssat[obs[j].sat-1].code[1];
            if (!test_sys(rtk->ssat[sati-1].sys,m,opt->qzsmodear,code)) continue;
            if (timediff(obs[0].time,nav->ssr[sati-1].t0[1])>=5.0) continue;
            l++;
        }
        for (f=opt->mode>PMODE_DGPS?0:nf;f<nf*2;f++) { /* freq*(phase+code) */
            /* search reference satellite with highest elevation */
            for (i=-1,j=0;j<n;j++) {
                flg=0;
                sati = obs[j].sat;
                sysi = rtk->ssat[sati-1].sys;
                code=rtk->ssat[obs[j].sat-1].code[f%nf];
                if (!test_sys(sysi,m,opt->qzsmodear,code)) continue;
                if (niter>0) {
                    for (k=0;k<niter;k++) {
                        if (refsat[NFREQ*2*MAXREF*m+NFREQ*2*k+f]==sati) flg=1;
                    }
                    if (flg==1) continue;
                }
                if (!validobs(j,f,nf,y)) continue;
                if (opt->qzsmodear==2&&opt->posopt[8]&&sysi==SYS_QZS) continue;
                if (l>1&&(timediff(obs[0].time,nav->ssr[sati-1].t0[1])>=5.0)) continue;
                if (i<0||azel[1+j*2]>=azel[1+i*2]) i=j;
            }
            if (i<0) continue;

            sati = obs[i].sat;
            refsat[NFREQ*2*MAXREF*m+NFREQ*2*niter+f]=sati;
            if (niter>0) {
                trace(2,"refsat changed %s L%d sat=%2d -> %2d\n",time_str(rtk->sol.time,0),
                      f%nf+1,refsat[NFREQ*MAXREF*2*m+NFREQ*2*(niter-1)+f],
                      refsat[NFREQ*MAXREF*2*m+NFREQ*2*niter+f]);
            }
            
            /* make single difference */
            for (j=0;j<n;j++) {
                if (i==j) continue;
                satj = obs[j].sat;
                sysi = rtk->ssat[sati-1].sys;
                sysj = rtk->ssat[satj-1].sys;
                code=rtk->ssat[obs[j].sat-1].code[f%nf];
                if (!test_sys(sysj,m,opt->qzsmodear,code)) continue;
                if (!validobs(j,f,nf,y)) continue;

                ff = f%nf;
                lami = nav->lam[sati-1][ff];
                lamj = nav->lam[satj-1][ff];
                if (lami<=0.0 || lamj<=0.0) continue;

                if (H) {
                    Hi=H+nv*rtk->nx;
                    for (k=0;k<rtk->nx;k++) Hi[k]=0.0;
                }

                v[nv]=y[f+i*nf*2]-y[f+j*nf*2]; /* single differenced residual */

                if (H) { /* partial derivatives by rover position */
                    for (k=0;k<3;k++) Hi[k]=-e[k+i*3]+e[k+j*3];
                }

                /* single-differenced ionospheric delay term */
                if (opt->ionoopt==IONOOPT_EST||opt->ionoopt==IONOOPT_EST_ADPT) {
                    fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
                    didxi=(f<nf?-1.0:1.0)*fi*fi*im[i];
                    didxj=(f<nf?-1.0:1.0)*fj*fj*im[j];
                    v[nv]-=didxi*x[II(sati,opt)]-didxj*x[II(satj,opt)];
                    if (H) {
                        Hi[II(sati,opt)]= didxi;
                        Hi[II(satj,opt)]=-didxj;
                    }
                }

                /* single-differenced tropospheric delay term */
                if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                    v[nv]-=(tropu[i]-tropu[j]);
                    for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {
                        if (!H) continue;
                        Hi[ITT(opt)+k]= (dtdxu[k+i*3]-dtdxu[k+j*3]);
                    }
                }
                /* single differenced phase-bias term */
                if (f<nf) {
                    v[nv] -= lami*x[IB(sati,f,opt)]-lamj*x[IB(satj,f,opt)];
                    if (H) {
                        Hi[IB(sati,f,opt)] = lami;
                        Hi[IB(satj,f,opt)] = -lamj;
                    }
                }

                /* test inovation*/
                /* deleted */

                if (f<nf) { /* phase */
                    rtk->ssat[satj-1].resc[f   ] = v[nv];
                    Ri[nv] = varerr(sati,sysi,azel[1+i*2],0,opt);
                    Rj[nv] = varerr(satj,sysj,azel[1+j*2],0,opt);
                    if (f==1) {
                        Ri[nv]*=SQR(2.55/1.55);
                        Rj[nv]*=SQR(2.55/1.55);
                    }
                } else {    /* code */
                    rtk->ssat[satj-1].resp[f-nf] = v[nv];
                    Ri[nv] = varerr(sati,sysi,azel[1+i*2],1,opt);
                    Rj[nv] = varerr(satj,sysj,azel[1+j*2],1,opt);
                }
                /* set valid data flags */
                if (f<nf) rtk->ssat[sati-1].vsat[f]=rtk->ssat[satj-1].vsat[f]=1;
                

                vflg[nv++] = (sati<<16)|(satj<<8)|((f<nf?0:1)<<4)|(f%nf);
                nb[b]++;
            }
            b++;
        }
    } /* end of system loop */
    
    /* measurement error covariance */
    ddcov(nb,b,Ri,Rj,nv,R);
    
    free(Ri);free(Rj);free(im);free(tropu);free(dtdxu);
    
    return nv;
}

/* single to double-difference transformation matrix (D') --------------------*/
static int ddmat(rtk_t *rtk, double *D, int *indI, int *indJ)
{
    int i,j,k,m,f,nb=0,nx=rtk->nx,na=rtk->na,nf=NF(&rtk->opt),code;

    trace(3,"ddmat   :\n");

    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].fix[j]=0;
    }
    for (i=0;i<na;i++) D[i+i*nx]=1.0;

    for (m=0;m<6;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds,4:qzs,5:gps(l2c) */
        if (m==1&&rtk->opt.glomodear==0) continue;
        if (m==4&&rtk->opt.qzsmodear!=1) continue;
        
        for (f=0,k=na;f<nf;f++,k+=MAXSAT) {
            for (i=k;i<k+MAXSAT;i++) {
                code=rtk->ssat[i-k].code[f];
                if (rtk->x[i]==0.0||
                    !test_sys(rtk->ssat[i-k].sys,m,rtk->opt.qzsmodear,code)||
                    !rtk->ssat[i-k].vsat[f]) {
                    continue;
                }
                if (rtk->ssat[i-k].lock[f]>0&&!(rtk->ssat[i-k].slip[f]&2)&&
                    rtk->ssat[i-k].azel[1]>=rtk->opt.elmaskar) {
                    rtk->ssat[i-k].fix[f]=2; /* fix */
                    break;
                }
                else rtk->ssat[i-k].fix[f]=1;
            }
            for (j=k;j<k+MAXSAT;j++) {
                code=rtk->ssat[j-k].code[f];
                if (i==j||rtk->x[j]==0.0||
                    !test_sys(rtk->ssat[j-k].sys,m,rtk->opt.qzsmodear,code)||
                    !rtk->ssat[j-k].vsat[f]) {
                    continue;
                }
                if (rtk->ssat[j-k].lock[f]>0&&!(rtk->ssat[j-k].slip[f]&2)&&
                    rtk->ssat[j-k].azel[1]>=rtk->opt.elmaskar) {

                    D[i+(na+nb)*nx]= 1.0;
                    D[j+(na+nb)*nx]=-1.0;
                    indI[nb]=i-na;indJ[nb]=j-na;
                    nb++;
                    rtk->ssat[j-k].fix[f]=2; /* fix */
                }
                else rtk->ssat[j-k].fix[f]=1;
            }
        }
    }

    trace(5,"D=\n"); tracemat(5,D,nx,na+nb,2,0);

    return nb;
}
/* restore single-differenced ambiguity --------------------------------------*/
static void restamb(rtk_t *rtk, const double *bias, int nb, double *xa)
{
    int i,n,m,f,index[MAXSAT],nv=0,nf=NF(&rtk->opt),code;
    
    trace(3,"restamb :\n");
    
    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i]; /* init all fixed states to float state values */
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i]; /* overwrite non phase-bias states with fixed values */
    
    for (m=0;m<6;m++) for (f=0;f<nf;f++) {
        
        for (n=i=0;i<MAXSAT;i++) {
            code=rtk->ssat[i].code[f];
            if (!test_sys(rtk->ssat[i].sys,m,rtk->opt.qzsmodear,code)
                ||rtk->ssat[i].fix[f]!=2) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
            rtk->ssat[i].lock[f]++;  /* increment count this sat used for AR */
        }
        if (n<2) continue;
        
        xa[index[0]]=rtk->x[index[0]];
        
        for (i=1;i<n;i++) {
            xa[index[i]]=xa[index[0]]-bias[nv++];
        }
    }
}
/* hold integer ambiguity ----------------------------------------------------*/
static void holdamb(rtk_t *rtk, const double *xa)
{
    double *v,*H,*R;
    int i,n,m,f,info,index[MAXSAT],nb=rtk->nx-rtk->na,nv=0,nf=NF(&rtk->opt),code;
    
    trace(3,"holdamb :\n");
    
    v=mat(nb,1); H=zeros(nb,rtk->nx);
    
    for (m=0;m<6;m++) for (f=0;f<nf;f++) {
        
        for (n=i=0;i<MAXSAT;i++) {
            code=rtk->ssat[i].code[f];
            if (!test_sys(rtk->ssat[i].sys,m,rtk->opt.qzsmodear,code)||
                rtk->ssat[i].fix[f]!=2||
                rtk->ssat[i].azel[1]<rtk->opt.elmaskhold||
                !rtk->ssat[i].vsat[f]) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
            rtk->ssat[i].fix[f]=3; /* hold */
        }
        /* constraint to fixed ambiguity */
        for (i=1;i<n;i++) {
            v[nv]=(xa[index[0]]-xa[index[i]])-(rtk->x[index[0]]-rtk->x[index[i]]);
            
            H[index[0]+nv*rtk->nx]= 1.0;
            H[index[i]+nv*rtk->nx]=-1.0;
            nv++;
        }
    }
    if (nv>0) {
        R=zeros(nv,nv);
        for (i=0;i<nv;i++) R[i+i*nv]=rtk->opt.varholdamb;
        
        /* update states with constraints */
        if ((info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv))) {
            trace(2,"filter error (info=%d)\n",info);
        }
        free(R);
    }
    free(v); free(H);
}

/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,ny,nb,nc,info,nx=rtk->nx,na=rtk->na,indI[MAXOBS*NFREQ],indJ[MAXOBS*NFREQ];
    double *D,*DP,*y,*b,*db,*Qb,*Qab,*QQ,s[2],thres;
    double *D_P;

    trace(3,"resamb_LAMBDA : nx=%d\n",nx);

    rtk->sol.ratio = 0.0;

    if (rtk->opt.mode <= PMODE_DGPS || rtk->opt.modear == ARMODE_OFF) {
        return 0;
    }

    /* single to double-difference transformation matrix (D') */
    D = zeros(nx,nx);
    nb=ddmat(rtk,D,indI,indJ);

    if (nb<=1) {
        trace(2,"no valid double-difference\n");
        free(D);
        return 0;
    }
    if (nb<rtk->opt.minamb) {
        trace(2,"less than minmum number of ambiguities nb=%2d\n",nb);
        free(D);
        return 0;
    }

    ny = na+nb; y=mat(ny,1); DP=mat(ny,nx);
    b=mat(nb,2); db=mat(nb,1); Qb=mat(nb,nb); Qab=mat(na,nb); QQ=mat(na,nb);
    nc=nx-na;D_P=mat(nb,nc);

    for (i=0;i<na;i++) {
        y[i]=rtk->x[i];
        for (j=0;j<nb;j++) {
            Qab[i+j*na]=rtk->P[i+(na+indI[j])*nx]-rtk->P[i+(na+indJ[j])*nx];
        }
    }
    for (i=0;i<nb;i++) {
        y[i+na]=rtk->x[na+indI[i]]-rtk->x[na+indJ[i]];
        for (j=0;j<nc;j++) {
            D_P[i+j*nb]=rtk->P[na+indI[i]+(na+j)*nx]-rtk->P[na+indJ[i]+(na+j)*nx]; 
        }
    }
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb[i+j*nb]=D_P[i+indI[j]*nb]-D_P[i+indJ[j]*nb];
    free(D_P);

    /* -----------------  conventional lambda   -----------------------  */
    if (!(info=lambda(nb,2,y+na,Qb,b,s))) {
        trace(5,"N(1)="); tracemat(5,b   ,1,nb,10,3);
        trace(5,"N(2)="); tracemat(5,b+nb,1,nb,10,3);

        rtk->sol.ratio = (s[0]>0) ? (float)(s[1]/s[0]) : 0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio = 999.9f;
        thres=qf[opt->alphaar][nb-1];

        /* validation by ratio-test */
        if (s[0]<=0.0||s[1]/s[0]>=thres) {
            for (i=0;i<na;i++) {
                rtk->xa[i] = rtk->x[i];
                for (j=0;j<na;j++) rtk->Pa[i+j*na] = rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i] = b[i];
                y[na+i] -= b[i];
            }
            if (!matinv(Qb,nb)) {
                /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
                matmul("NN",nb,1,nb, 1.0,Qb ,y+na,0.0,db);
                matmul("NN",na,1,nb,-1.0,Qab,db  ,1.0,rtk->xa);

                /* Covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
                matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,rtk->Pa);

                trace(3,"resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                      nb, (s[0]==0.0) ? 0.0 : s[1]/s[0], s[1], s[0]);

                /* restore single-differenced ambiguity */
                restamb(rtk,bias,nb,xa);
            } else {
                trace(2, "lambda error nb=%d\n", nb);
                nb = 0;
            }
        } else { /* validation failed */
            trace(3,"ambiguity validation failed (nb=%d thres=%.2f ratio=%.2f s=%.2f/%.2f)\n",
                  nb, thres, (s[0]==0.0) ? 0.0 : s[1]/s[0], s[1], s[0]);
            nb = 0;
        }
    } else {
        trace(2,"lambda error (info=%d)\n",info);
        nb = 0;
    }

    free(D); free(y); free(DP);
    free(b); free(db); free(Qb); free(Qab); free(QQ);

    return nb; /* number of ambiguities */
}

static void elsort(rtk_t *rtk, const obsd_t *obs, int *isat, double *el,
                   const int n)
{
    double el_;
    int i,j,k,m=0,sati;

    for (i=0;i<n;i++) {
        sati=obs[i].sat;
        el_=rtk->ssat[sati-1].azel[1];
        if (m<=0) {
            isat[m] = sati;
            el[m++] = el_;
        } else {
            for (j=0;j<m;j++) if (el_<el[j]) break;
            if (j>=n) continue;
            for (k=m; k>j;k--) {
                isat[k] = isat[k-1];
                el[k] = el[k-1];
            }
            isat[j] = sati;
            el[j] = el_;
            if (m<n) m++;
        }
    }
}

/* compute single difference dop */
static double sddop(const rtk_t *rtk, const double *H, const int *vflg, 
                    const double elmask, const int nv)
{
    int i,cnt,freq,type,sat;
    double *H_,Q[9],dop=1000.0;
    
    H_=mat(3,nv);

    for (i=cnt=0;i<nv;i++) {
        sat=(vflg[i]>>8)  & 0xFF;
        type=(vflg[i]>>4) & 0xF;
        freq=(vflg[i])    & 0xF;
        if (freq>0||type>0||rtk->ssat[sat-1].vsat[0]==0||
                    rtk->ssat[sat-1].azel[1]<=elmask) continue;
        H_[cnt*3  ]=H[i*rtk->nx];
        H_[cnt*3+1]=H[i*rtk->nx+1];
        H_[cnt*3+2]=H[i*rtk->nx+2];
        cnt++;
    }
    matmul("NT",3,3,cnt,1.0,H_,H_,0.0,Q);
    if (cnt>0&&!matinv(Q,3)) {
        dop=sqrt(Q[0]+Q[4]+Q[8]); /*PDOP*/
    }
    
    free(H_);
    
    return dop;
}

/* compute zero difference (conventional) dop */
static double zddop(const rtk_t *rtk, const double *azel,
                    const obsd_t *obs, const double elmask, const int n)
{
    int i,j;
    double *azel_,dop[4];
    
    azel_=mat(2,n);
    
    for (i=j=0;i<n;i++) {
        if (rtk->ssat[obs[i].sat-1].vsat[0]==0) continue;
        azel_[j*2  ]=azel[i*2  ];
        azel_[j*2+1]=azel[i*2+1];
        j++;
    }
    
    dops(j,azel_,elmask,dop);
    
    free(azel_);
    
    return dop[1];  /* PDOP */
}

/* compute dop */
static double calpdop(const rtk_t *rtk, const double *H, const double *azel,
                      const int *vflg, const obsd_t *obs, const double elmask, 
                      const int n, const int nv)
{
    switch (rtk->opt.refdop) {
        case 0: return zddop(rtk,azel,obs,elmask,n); break;
        case 1: return sddop(rtk,H,vflg,elmask,nv);  break;
    }
    
    return 0.0;
}

/* number of estimated states ------------------------------------------------*/
extern int ppp_rtk_nx(const prcopt_t *opt)
{
    return NX(opt);
}

static void resetfilterflag(nav_t *nav)
{
    nav->ionoreset = FALSE;
    nav->filreset = FALSE;
    clearsatcorr();
}
/* precise positioning ---------------------------------------------------------
* input observation data and navigation message, compute rover position by
* precise positioning
* args   : rtk_t *rtk       IO  rtk control/result struct
*            rtk->sol       IO  solution
*                .time      O   solution time
*                .rr[]      IO  rover position/velocity
*                               (I:fixed mode,O:single mode)
*                .dtr[0]    O   receiver clock bias (s)
*                .dtr[1]    O   receiver glonass-gps time offset (s)
*                .Qr[]      O   rover position covarinace
*                .stat      O   solution status (SOLQ_???)
*                .ns        O   number of valid satellites
*                .age       O   age of differential (s)
*                .ratio     O   ratio factor for ambiguity validation
*            rtk->rb[]      IO  base station position/velocity
*                               (I:relative mode,O:moving-base mode)
*            rtk->nx        I   number of all states
*            rtk->na        I   number of integer states
*            rtk->ns        O   number of valid satellite
*            rtk->tt        O   time difference between current and previous (s)
*            rtk->x[]       IO  float states pre-filter and post-filter
*            rtk->P[]       IO  float covariance pre-filter and post-filter
*            rtk->xa[]      O   fixed states after AR
*            rtk->Pa[]      O   fixed covariance after AR
*            rtk->ssat[s]   IO  sat(s+1) status
*                .sys       O   system (SYS_???)
*                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)
*                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)
*                .vs   [r]  O   data valid single     (r=0:rover,1:base)
*                .resp [f]  O   freq(f+1) pseudorange residual (m)
*                .resc [f]  O   freq(f+1) carrier-phase residual (m)
*                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)
*                .fix  [f]  O   freq(f+1) ambiguity flag
*                               (0:nodata,1:float,2:fix,3:hold)
*                .slip [f]  O   freq(f+1) slip flag
*                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
*                                bit2:parity unknown, bit1:slip)
*                .lock [f]  IO  freq(f+1) carrier lock count
*                .outc [f]  IO  freq(f+1) carrier outage count
*                .slipc[f]  IO  freq(f+1) cycle slip count
*                .rejc [f]  IO  freq(f+1) data reject count
*                .gf        IO  geometry-free phase (L1-L2) (m)
*                .gf2       IO  geometry-free phase (L1-L5) (m)
*            rtk->nfix      IO  number of continuous fixes of ambiguity
*            rtk->neb       IO  bytes of error message buffer
*            rtk->errbuf    IO  error message buffer
*            rtk->tstr      O   time string for debug
*            rtk->opt       I   processing options
*          obsd_t *obs      I   observation data for an epoch
*                               obs[i].rcv=1:rover,2:reference
*                               sorted by receiver and satellte
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation messages
* return : status (0:no solution,1:valid solution)
* notes  : before calling function, base station position rtk->sol.rb[] should
*          be properly set for relative mode except for moving-baseline
*-----------------------------------------------------------------------------*/
#ifndef CSSR2OSR_VRS
extern void ppp_rtk_pos(rtk_t *rtk, const obsd_t *obs, int n, nav_t *nav)
{
    static gtime_t regularly = {-1, 0.0};
    prcopt_t *opt=&rtk->opt;
    static grid_t grid;
    double *rs,*dts,*var,*y,*e,*v,*azel,*H,*R,*xp,*Pp,*Qp,*xa,*bias;
    double el[MAXSAT];
    int i,j,k=0,l,f,nv,info,nf=rtk->opt.nf,sati;
    int stat = (rtk->opt.mode <= PMODE_DGPS) ? SOLQ_DGPS : SOLQ_FLOAT;
    int vflg[MAXOBS*NFREQ*4+1],svh[MAXOBS],nb,isat[MAXSAT];
    float age;
    static double cpc[MAXSAT*NFREQ]={0};
    static gtime_t pt0[MAXSAT]={0};
    double cpc_[MAXSAT*NFREQ],dist;
    gtime_t pt0_[MAXSAT];
    static int resetcnt=0,retrycnt=-1000,cntdiffp=0,np;
    double pdop;
    static int backup = FALSE;
    double pos[3];
    gtime_t temp;
    int nn;

    trace(2,"ppp_rtk_pos   : time=%s nx=%d n=%d filreset=%d prevstat=%d\n",time_str(obs[0].time,0),rtk->nx,n,nav->filreset,rtk->sol.pstat);

    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);

    nv=n*nf*2*3;

    y = mat(nv,1); e = mat(n,3);

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys = satsys(i+1,NULL);
        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat[j] = rtk->ssat[i].snr[j] = 0;
    }

    for (i=0;i<MAXSAT;i++) rtk->ssat[i].fix[0]=0;

    /* reset all states if time difference between current
    and previous exceeds MAXOBSLOSS */
    if (rtk->tt>opt->maxobsloss) {
        for (i=0;i<rtk->nx;i++) rtk->x[i]=0.0;
    }
    
    /* temporal update of states */
    udstate_ppp(rtk,obs,n,nav);

    grid.network = opt->netnum;
    ecef2pos(rtk->x, pos);
    if ((nn = get_grid_index(nav, pos, &grid, opt, obs[0].time, FALSE)) == 0) {
        if (backup == FALSE || timediff(obs[0].time, (temp = GetBackupCSSRTime())) > 180.0) {
            free(rs); free(dts); free(var); free(y); free(e); free(azel);
            rtk->sol.stat = SOLQ_SINGLE;
            return;
        } else {
            trace(3, "use backup CSSR data: obs=%.1f, CSSR=%.1f, age=%.1f\n",
                time2gpst(obs[0].time, NULL), time2gpst(temp, NULL),
                timediff(obs[0].time, temp));
            RestoreCurrentCSSR(&grid);
        }
    } else {
        if (GetCloseCSSR(obs[0].time, grid.network) == FALSE) {
            free(rs); free(dts); free(var); free(y); free(e); free(azel);
            rtk->sol.stat = SOLQ_SINGLE;
            return;
        }
        BackupCurrentCSSR(&grid);
        backup = TRUE;
    }
    
    nav->facility = GetCurrentCSSRFacility();
    nav->ssrtime = GetCurrentCSSRTime();
    CheckCSSRFacility(nav, grid.network);
    
    if (opt->posopt[11] != 0) {
        if (!opt->posopt[9] && IsSISAdjust()) {
            rtk->sisadjust = 1;
        } else {
            rtk->sisadjust = 0;
        }
    } else {
        rtk->sisadjust = IsSISAdjust();
    }
    for (i = 0; i < MAXSAT; ++i) {
        UpdateGlobalCSSR(&nav->ssr[i], i + 1);
    }
    UpdateLocalCSSR(nav);

    
    regularly = (regularly.time == -1 ? obs[0].time: regularly);
    if (opt->regularly != 0 && timediff(obs[0].time, regularly) >= (double)opt->regularly) {
        double tow = time2gpst(timeget(), NULL);
        trace(1, "ppp_rtk_pos(): reset filter, tow=%.2f, network=%d\n", tow, rtk->opt.netnum);
        for (i=NP(opt);i<rtk->nx;i++) rtk->x[i]=0.0;
        regularly = obs[0].time;
        resetfilterflag(nav);
        stat = SOLQ_NONE;
        rtk->sol.stat=SOLQ_SINGLE;
        ClearCurrentCSSR();
        retrycnt=opt->retrycnt;
        resetcnt=0;
        udstate_ppp(rtk,obs,n,nav);
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return;
    }
    if (nav->filreset == TRUE ||
        (opt->epochtoretry>0&&retrycnt>0&&rtk->sol.pstat!=1&&resetcnt==opt->epochtoretry)) {
        double tow = time2gpst(timeget(), NULL);
        trace(1, "ppp_rtk_pos(): reset filter, tow=%.2f, network=%d, filreset=%d\n", tow, rtk->opt.netnum, nav->filreset);
        for (i=NP(opt);i<rtk->nx;i++) rtk->x[i]=0.0;
        resetfilterflag(nav);
        stat = SOLQ_FLOAT;
        rtk->sol.stat=SOLQ_FLOAT;
        if (nav->filreset == TRUE) {
            retrycnt=opt->retrycnt;
        } else {
            --retrycnt;
        }
        resetcnt=0;
        udstate_ppp(rtk,obs,n,nav);
    }
    
    if (nav->ionoreset == TRUE) {
        trace(2, "change select network or grid: tow=%.1f\n", time2gpst(obs[0].time, NULL));
        for (i = 1; i <= MAXSAT; ++i) {
            rtk->x[II(i,opt)] = 0.0;
        }
        udstate_ppp(rtk, obs, n, nav);
        nav->ionoreset = FALSE;
    }
    ++resetcnt;
    trace(4,"x(0)="); tracemat(4,rtk->x,1,NR(opt),13,4);

    /* satellite positions and clocks */
    saveposition(rtk->x);
    satposs(obs[0].time,obs,n,nav,rtk->opt.sateph,rs,dts,var,svh);
    for (i=0;i<MAXSAT;i++) {
        satpos_ssr_sis(obs[0].time,obs[0].time,i+1,rtk,nav);
    }

    xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx); xa = mat(rtk->nx,1);
    Qp=zeros(rtk->nx,rtk->nx);
    
    matcpy(xp,rtk->x,rtk->nx,1);

    v=mat(nv,1);  R=mat(nv,nv); H=mat(rtk->nx,nv);
    bias = mat(rtk->nx,1);

    for (i=0;i<rtk->opt.niter;i++) {
        for (k=0;k<MAXREF;k++) {
            for (l=0;l<MAXSAT;l++) {
                pt0_[l]=pt0[l];
                for (f=0;f<NFREQ;f++) cpc_[l*NFREQ+f]=cpc[l*NFREQ+f];
            }
            nv=zdres(obs,n,rs,dts,var,svh,nav,xp,y,e,azel,rtk,TRUE,cpc_,pt0_,&grid,rtk->ssat,&rtk->opt,&rtk->sol,NULL);
            /* trop estimation params is set to be zero when trop ssr is valid */
            if ((resetcnt==1)&&rtk->opt.tropopt>=TROPOPT_EST&&nav->invtrop==0) {
                j=ITT(&rtk->opt);
                initx(rtk,INIT_ZWD,0.0,j);
                rtk->opt.std[2]=rtk->opt.prn[2]=0.0;
                if (rtk->opt.tropopt==TROPOPT_ESTG) {
                    for (l=0;l<2;l++) initx(rtk,1E-6,0.0,++j);
                }
            }
            if (!nv) {
                trace(2,"rover initial position error\n");
                stat=SOLQ_NONE;
                break;
            }
            /* double-differenced residuals and partial derivatives */
            if ((nv = ddres(rtk,nav,xp,Pp,obs,y,e,azel,n,v,H,R,vflg,k))<=0) {
                trace(2,"no double-differenced residual\n");
                stat=SOLQ_NONE;
                break;
            }
            /* Kalman filter measurement update */
            matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
            if ((info=filter2(rtk,xp,Pp,Qp,H,v,R,rtk->nx,nv,vflg,0))!=0) {
                trace(2,"ppp-rtk filter error %s info=%d\n",time_str(rtk->sol.time,0),info);
                continue;
            }
            /* postfit residuals of float solution */
            if (zdres(obs,n,rs,dts,var,svh,nav,xp,y,e,azel,rtk,FALSE,cpc_,pt0_,&grid,rtk->ssat,&rtk->opt,&rtk->sol,NULL)) {
                nv = ddres(rtk,nav,xp,Pp,obs,y,e,azel,n,v,NULL,R,vflg,k);
                /* validation of float solution */
                filter2(rtk,xp,Pp,NULL,H,v,R,rtk->nx,nv,vflg,1);
                if (rtk->sol.chisq<100.0) break;
            }
            matcpy(xp,rtk->x,rtk->nx,1);
            if (k==(MAXREF-1)) {
                stat=SOLQ_NONE;
                trace(2,"measurement update failed \n");
            }
        }
    }

    /* ssr age */
    rtk->sol.age=1e4;
    for (i=0;i<n&&i<MAXOBS;i++) {
        sati = obs[i].sat;
        age=timediff(obs[i].time,nav->ssr[sati-1].t0[1]);
        if (rtk->sol.age>age) rtk->sol.age=age;
    }

    /* update process noise of ionoshepre delay */
    if (rtk->opt.ionoopt==IONOOPT_EST_ADPT) {
        for (i=0;i<n;i++) {
            j=II(obs[i].sat,&rtk->opt);
            rtk->Q[j+j*rtk->nx]=rtk->opt.forgetion*rtk->Q[j+j*rtk->nx]+
                (1.0-rtk->opt.forgetion)*SQR(rtk->opt.afgainion)*Qp[j+j*rtk->nx];
        }
    }
    /* update process noise of pos, vel, and acc */
    if (rtk->opt.prnadpt) {
        np=NP(opt);
        for (i=0;i<np;i++) for (j=0;j<np;j++) {
            rtk->Q[i+j*rtk->nx]=rtk->opt.forgetpva*rtk->Q[i+j*rtk->nx]+
                (1.0-rtk->opt.forgetpva)*SQR(rtk->opt.afgainpva)*Qp[i+j*rtk->nx];
        }
    }
    
    if (stat!=SOLQ_NONE) {
        /* upate state and covariance matrix */
        matcpy(rtk->x, xp, rtk->nx, 1);
        matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
        /* update ambiguity control struct */
        rtk->sol.ns = 0;
        for (i=0;i<n;i++) {
            sati = obs[i].sat;
            for (f=0;f<nf;f++) {
                if (!rtk->ssat[sati-1].vsat[f]) continue;
                rtk->ssat[sati-1].lock[f]++;
                rtk->ssat[sati-1].outc[f]=0;
                if (f==0) rtk->sol.ns++; /* valid satellite count by L1 */
            }
        }
        /* lack of valid satellites */
        if (rtk->sol.ns<4) {
            trace(2,"lack of valid satellites.\n");
            stat=SOLQ_NONE;
        }
    }

    /* keep symmetry of covariance matrix */
    for (i=0;i<rtk->nx;i++) for (j=0;j<i;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=(rtk->P[i+j*rtk->nx]+rtk->P[j+i*rtk->nx])*0.5;
    }

    /* resolve integer ambiguity by LAMBDA */
    pdop=calpdop(rtk,H,azel,vflg,obs,rtk->opt.elmaskar,n,nv);
    if (stat!=SOLQ_NONE&&(rtk->opt.maxpdopar==0.0||pdop<=rtk->opt.maxpdopar)) {
        if ((nb=resamb_LAMBDA(rtk,bias,xa))<=1&&rtk->opt.posopt[6]) {
            int exsat=0;
            double ratio=0.0;
            elsort(rtk,obs,isat,el,n);
            for (l=0;l<rtk->opt.armaxdelsat;l++) {
                for (i=0;i<n;i++) {
                    sati = isat[i];
                    if (el[i]<=rtk->opt.elmaskar) continue;
                    for (f=j=0;f<nf;f++) {
                        if (!rtk->ssat[sati-1].vsat[f]||rtk->ssat[sati-1].lock[f]<=0) j++;
                    }
                    if (j>1) continue;
                    for (f=0;f<nf;f++) rtk->ssat[sati-1].vsat[f]=0;
                    pdop=calpdop(rtk,H,azel,vflg,obs,rtk->opt.elmaskar,n,nv);
                    trace(2, "exclude sat: sys=%d, sat=%d, l=%d                    \n", satsys(sati,NULL), sati, l);
                    if ((rtk->opt.maxpdopar==0.0||pdop<=rtk->opt.maxpdopar)&&
                                        (nb=resamb_LAMBDA(rtk,bias,xa))>1) {
                        trace(2,"PAR OK l=%d i=%2d nb=%2d obs=%.1f excluded sat=%2d el=%4.1f pdop=%4.1f ratio=%5.2f\n",
                              l,i,nb,time2gpst(obs[0].time,NULL),sati,el[i]*R2D,pdop,rtk->sol.ratio);

                        break;
                    } else {
                        for (f=0;f<nf;f++) rtk->ssat[sati-1].vsat[f]=1;
                        if (ratio<rtk->sol.ratio) {
                            exsat=sati;ratio=rtk->sol.ratio;
                        }
                    }
                }
                if (nb>1) {
                    break;
                } else {
                    for (f=0;f<nf;f++) rtk->ssat[exsat-1].vsat[f]=0;
                }
            }
        }
        
        if (nb>1&&zdres(obs,n,rs,dts,var,svh,nav,xa,y,e,azel,rtk,FALSE,cpc_,pt0_,&grid,rtk->ssat,&rtk->opt,&rtk->sol,NULL)) {
            /* post-fix residuals for fixed solution */
            nv = ddres(rtk,nav,xa,NULL,obs,y,e,azel,n,v,NULL,R,vflg,k);
            /* validation of fixed solution */
            filter2(rtk,xp,Pp,NULL,H,v,R,rtk->nx,nv,vflg,2);
            
            if (rtk->sol.chisq<opt->maxinno[3]) {
                /* hold integer ambiguity */
                pdop=calpdop(rtk,H,azel,vflg,obs,rtk->opt.elmaskhold,n,nv);
                if (++rtk->nfix >= rtk->opt.minfix &&
                    rtk->opt.modear == ARMODE_FIXHOLD &&
                    (rtk->opt.maxpdophold==0.0||pdop<=rtk->opt.maxpdophold)) {
                    holdamb(rtk,xa);
                }
           }
           if (rtk->sol.chisq<opt->maxinno[4]) stat = SOLQ_FIX;
        }
    }
    
    /* reset all states due to large differece in solution between point pos and pprtk solutions */
    if (rtk->opt.maxdiffp>0.0&&norm(rtk->x,3)>0.0&&norm(rtk->sol.rr,3)>0.0) {
        dist=SQR(rtk->x[0]-rtk->sol.rr[0])+SQR(rtk->x[1]-rtk->sol.rr[1])+SQR(rtk->x[2]-rtk->sol.rr[2]);
        if (dist>SQR(rtk->opt.maxdiffp)) {
            ++cntdiffp;
            if (cntdiffp>rtk->opt.poserrcnt) {
                trace(2,"reset all states due to large pos error dist = %.2f\n",sqrt(dist));
                stat=SOLQ_NONE;
                for (i=0;i<3;i++) {
                    rtk->x[i]=rtk->sol.rr[i] ;
                    rtk->P[i+i*rtk->nx]=(double)rtk->sol.qr[i];
                }
                cntdiffp=0;
            }
        } else {
            cntdiffp=0;
        }
    }
    
    /* save solution status */
    if (stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i] = rtk->xa[i];
            rtk->sol.qr[i] = (float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3] = (float)rtk->Pa[1];
        rtk->sol.qr[4] = (float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5] = (float)rtk->Pa[2];
    } else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i] = rtk->x[i];
            rtk->sol.qr[i] = (float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3] = (float)rtk->P[1];
        rtk->sol.qr[4] = (float)rtk->P[1+2*rtk->nx];
        rtk->sol.qr[5] = (float)rtk->P[2];
        rtk->nfix = 0;
    }

    for (i=0; i<n; i++) {
        for (j=0; j<nf; j++) {
            if (obs[i].L[j]==0.0) continue;
            rtk->ssat[obs[i].sat-1].pt[obs[i].rcv-1][j] = obs[i].time;
            rtk->ssat[obs[i].sat-1].ph[obs[i].rcv-1][j] = obs[i].L[j];
        }
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nf;j++) {
            rtk->ssat[obs[i].sat-1].snr[j] = obs[i].SNR[j];
        }
    }
    for (i=0;i<MAXSAT;i++) {
        for (j=0;j<nf;j++) {
            if (rtk->ssat[i].fix[j]==2 && stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
            if (rtk->ssat[i].slip[j]&1) rtk->ssat[i].slipc[j]=1;
        }
    }

    if (stat!=SOLQ_NONE) rtk->sol.stat=stat;
    else rtk->sol.stat=SOLQ_SINGLE;

    if (stat==SOLQ_NONE) {
        for (i=0;i<rtk->nx;i++) rtk->x[i]=0.0;
    }
    
    for (l=0;l<MAXSAT;l++) {
        pt0[l]=pt0_[l];
        for (f=0;f<NFREQ;f++) cpc[l*NFREQ+f]=cpc_[l*NFREQ+f];
    }
    
    rtk->sol.pstat=stat;
    dops(n,azel,opt->elmin,rtk->sol.dop); /* {GDOP,PDOP,HDOP,VDOP} */

    free(rs); free(dts); free(var); free(y) ; free(e); free(azel);
    free(xp); free(Pp); free(xa); free(v); free(H); free(R); free(bias);
    free(Qp);
}
#endif /* CSSR2OSR_VRS */
