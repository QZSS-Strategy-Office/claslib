/*------------------------------------------------------------------------------
* rtkvrs.c : post-processing positioning for VRS-RTK
*            Derived from postpos.c
*
*          Copyright (C) 2007- by T.TAKASU, All rights reserved.
*          Copyright (C) 2015- by Mitsubishi Electric Corp., All rights reserved.
*
* version : $Revision:  1.0
* history   2019/2/ 1   1.0 new
*
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/
#define MAXREF 2                    /* max interation of searching reference sat */

#define SQR(x)      ((x)*(x))
#define MAXITER 20            /* max iteration of PAR */
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS     SQR(30.0)       /* initial variance of receiver pos ((m)^2) */
#define VAR_VEL     SQR(1.0)        /* initial variance of receiver acc ((m/s)^2) */
#define VAR_ACC     SQR(1.0)        /* initial variance of receiver acc ((m/ss)^2) */
#define VAR_HWBIAS  SQR(1.0)  /* initial variance of h/w bias ((m/MHz)^2) */
#define VAR_GRA     SQR(0.001) /* initial variance of gradient (m^2) */
#define INIT_ZWD    0.15     /* initial zwd (m) */

#define PRN_HWBIAS  1E-6     /* process noise of h/w bias (m/MHz/sqrt(s)) */
#define GAP_RESION  120      /* gap to reset ionosphere parameters (epochs) */
#define MAXACC      30.0     /* max accel for doppler slip detection (m/s^2) */

#define VAR_HOLDAMB 0.001    /* constraint to hold ambiguity (cycle^2) */

#define TTOL_MOVEB  (1.0+2*DTTOL)
                             /* time sync tolerance for moving-baseline (s) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((((opt)->ionoopt!=IONOOPT_EST)&&((opt)->ionoopt!=IONOOPT_EST_ADPT))?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */
#define ITT(opt)    (NP(opt)+NI(opt))  /* tropos (r:0=rov,1:ref) */

enum ddres_order {
    FST = 1,
    SND,
    TRD
};

/* global variables ----------------------------------------------------------*/

static const double qfrtk[6][60] = {
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

/* save error message --------------------------------------------------------*/
static void errmsg(rtk_t *rtk, const char *format, ...)
{
    char buff[256],tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time,tstr,2);
    n=sprintf(buff,"%s: ",tstr+11);
    va_start(ap,format);
    n+=vsprintf(buff+n,format,ap);
    va_end(ap);
    n=n<MAXERRMSG-rtk->neb?n:MAXERRMSG-rtk->neb;
    memcpy(rtk->errbuf+rtk->neb,buff,n);
    rtk->neb+=n;
    trace(2,"%s",buff);
}

/* single-differenced measurement error variance -----------------------------*/
static double varerr(int sat, int sys, double el, double bl, double dt, int f,
                     const prcopt_t *opt)
{
    double a,b,c=opt->err[3]*bl/1E4,d=CLIGHT*opt->sclkstab*dt,fact=1.0;
    double sinel=sin(el);
    int i=sys==SYS_GLO?1:(sys==SYS_GAL?2:0),nf=NF(opt);
    
    /* extended error model */
    if (f>=nf&&opt->exterr.ena[0]) { /* code */
        a=opt->exterr.cerr[i][  (f-nf)*2];
        b=opt->exterr.cerr[i][1+(f-nf)*2];
        if (sys==SYS_SBS) {a*=EFACT_SBS; b*=EFACT_SBS;}
    }
    else if (f<nf&&opt->exterr.ena[1]) { /* phase */
        a=opt->exterr.perr[i][  f*2];
        b=opt->exterr.perr[i][1+f*2];
        if (sys==SYS_SBS) {a*=EFACT_SBS; b*=EFACT_SBS;}
    }
    else { /* normal error model */
        if (f>=nf) fact=opt->eratio[f-nf];
        if (fact<=0.0)  fact=opt->eratio[0];
        fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
        a=fact*opt->err[1];
        b=fact*opt->err[2];
        c=d=0.0;
    }

   return 1.0*(opt->ionoopt==IONOOPT_IFLC?3.0:1.0)*(a*a+b*b/sinel/sinel+c*c)+d*d;

}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
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
/* select common satellites between rover and reference station --------------*/
static int selsat(const obsd_t *obs, double *azel, int nu, int *nr,
                  const prcopt_t *opt, int *sat, int *iu, int *ir, int ch)
{
    int i,j,k=0,st,end;
    
    trace(3,"selsat  : nu=%d nr=%d\n",nu,nr);
    for(i=0;i<MAXSAT;i++) { sat[i]=-1;ir[i]=-1;iu[i]=-1;}

    st=(ch==0)?nu:nu+nr[0];
    end=(ch==0)?nu+nr[0]:nu+nr[0]+nr[1];    
    for (i=0,j=st;i<nu&&j<end;i++,j++) {
        if      (obs[i].sat<obs[j].sat) j--;
        else if (obs[i].sat>obs[j].sat) i--;
        else if (azel[1+j*2]>=opt->elmin) { /* elevation at base station */
            if (timediff(obs[i].time,obs[j].time)<-1.0) {
                trace(5,"ch:%d obs[i].sat:%d dt:%f\n",ch,obs[i].sat,timediff(obs[i].time,obs[j].time));
                continue;
            }
            sat[k]=obs[i].sat; iu[k]=i; ir[k++]=j;
            trace(4,"(%2d) sat=%3d iu=%2d ir=%2d\n",k-1,obs[i].sat,i,j);
        }
    }
    return k;
}
static int selfreqpair(const int sat, const prcopt_t *opt,const obsd_t *obs)
{
    int optf=opt->posopt[10];
    if (NFREQ==1||optf==POSL1) return 0;
    if (NFREQ>=3) {
        if ((satsys(sat,NULL))&(SYS_GAL)) return 2; /* GAL */
        if (optf==POSL1L2L5) return 1+2;
        if (optf==POSL1L5) return 2;
        if (optf==POSL1L5_L2&&obs->L[2]!=0.0&&obs->P[2]!=0.0) return 2;
    }
    return 1;
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void udpos(rtk_t *rtk, const obsd_t *obs, double tt)
{
    double *F,*FP,*xp,pos[3],Q[9]={0},Qv[9],var=0.0;
    int i, j, na,nb,nx,ni,flag_init=0;
    double *Fx,*Pxx,*Pxi,*Pxb,*Pxi_,*Pxb_,*xp_,ki;

    
    trace(3,"udpos   : tt=%.3f\n",tt);
    
    nx=rtk->nx;

    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
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
    for (i=0;i<3;i++) { var+=rtk->P[i+i*nx]; } var/=3.0;
    
    if (var>VAR_POS) {
        /* reset position with large variance */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    
    /* state transition of position/velocity/acceleration */
    if (rtk->opt.ionoopt>=IONOOPT_EST){
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
        F=eye(9); xp=mat(9,1); FP=mat(nx,nx);
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
/* temporal update of ionospheric parameters ---------------------------------*/
static void udion(rtk_t *rtk, double tt, double bl, const obsd_t *obs, int ns)
{
    double el,fact,qi;
    int i,j,f,k,m;

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
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop(rtk_t *rtk, double tt, double bl)
{
    int i,j,k;
    
    trace(3,"udtrop  : tt=%.1f\n",tt);
    
    for (i=0;i<2;i++) {
        j=IT(i,&rtk->opt);
        
        if (rtk->x[j]==0.0) {
            initx(rtk,INIT_ZWD,SQR(rtk->opt.std[2]),j); /* initial zwd */
            
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
/* temporal update of receiver h/w biases ------------------------------------*/
static void udrcvbias(rtk_t *rtk, double tt)
{
    int i,j;
    
    trace(3,"udrcvbias: tt=%.1f\n",tt);
    
    for (i=0;i<NFREQGLO;i++) {
        j=IL(i,&rtk->opt);
        
        if (rtk->x[j]==0.0) {
            initx(rtk,1E-6,VAR_HWBIAS,j);
        }
        /* hold to fixed solution */
        else if (rtk->nfix>=rtk->opt.minfix&&rtk->sol.ratio>rtk->opt.thresar[0]) {
            initx(rtk,rtk->xa[j],rtk->Pa[j+j*rtk->na],j);
        }
        else {
            rtk->P[j+j*rtk->nx]+=SQR(PRN_HWBIAS)*tt;
        }
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
static void udstate(rtk_t *rtk, const obsd_t *obs, int nu, const nav_t *nav)
{
    double tt=fabs(rtk->tt),bl=0.0,dr[3];
    
    trace(3,"udstate : nu=%d\n",nu);
    
    /* temporal update of position/velocity/acceleration */
    udpos(rtk,obs,tt);
    
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt>=TROPOPT_EST) {
        udtrop(rtk,tt,bl);
    }
    /* temporal update of phase-bias */
    if (rtk->opt.mode>PMODE_DGPS) {
        udbias_ppp(rtk,tt,obs,nu,nav);
    }
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        bl=baseline(rtk->x,rtk->rb,dr);
        udion(rtk,tt,bl,obs,nu);
    }
    /* temporal update of eceiver h/w bias */
    if (rtk->opt.glomodear==2&&(rtk->opt.navsys&SYS_GLO)) {
        udrcvbias(rtk,tt);
    }

}
/* undifferenced phase/code residual for satellite ---------------------------*/
static void zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
                      const double *azel, const double *dant,
                      const prcopt_t *opt, double *y)
{
    const double *lam=nav->lam[obs->sat-1];
    int i,nf=NF(opt),qj;
    
    qj=selfreqpair(obs->sat,opt,obs);
    for (i=0;i<nf;i++) {
        if (lam[i]==0.0) continue;
        
        /* check snr mask */
        if (testsnr(base,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) {
            continue;
        }
        if (i>0&&!(i&qj)) continue;
        /* residuals = observable - pseudorange */
        if (obs->L[i]!=0.0) y[i   ]=obs->L[i]*lam[i]-r-dant[i];
        if (obs->P[i]!=0.0) y[i+nf]=obs->P[i]       -r-dant[i];
    }
    
}
/* undifferenced phase/code residuals ----------------------------------------*/
static int zdres(rtk_t *rtk, int base, const obsd_t *obs, int n, const double *rs,
                 const double *dts, const int *svh, const nav_t *nav,
                 const double *rr, const prcopt_t *opt, int index, double *y,
                 double *e, double *azel)
{
    double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
    double zhd,zazel[]={0.0,90.0*D2R},ydif=0.0;
    int i,j,nf=NF(opt);
    
    trace(3,"zdres   : n=%d\n",n);
    
    for (i=0;i<n*nf*2;i++) y[i]=0.0;
    
    if (norm(rr,3)<=0.0) return 0; /* no receiver position */
    
    for (i=0;i<3;i++) rr_[i]=rr[i];
    
    /* earth tide correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rr_,opt->tidecorr,&nav->erp,
                 opt->odisp[base],NULL,0,0,NULL,NULL,NULL,NULL,disp);
        for (i=0;i<3;i++) rr_[i]+=disp[i];
    }
    ecef2pos(rr_,pos);
    
    for (i=0;i<n;i++) {
        /* compute geometric-range and azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr_,e+i*3))<=0.0) continue;
        if (satazel(pos,e+i*3,azel+i*2)<opt->elmin) continue;
        
        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;
        
        /* satellite clock-bias */
        r+=-CLIGHT*dts[i*2];
        
        /* troposphere delay model (hydrostatic) */
        zhd=tropmodel(obs[0].time,pos,zazel,0.0);
        r+=tropmapf(obs[i].time,pos,azel+i*2,NULL)*zhd;
        
        /* receiver antenna phase center correction */
        antmodel(opt->pcvr+index,opt->antdel[index],azel+i*2,opt->posopt[1],
                 dant);
        
        /* undifferenced phase/code residual for satellite */
        zdres_sat(base,r,obs+i,nav,azel+i*2,dant,opt,y+i*nf*2);
        
        /* 1/4 cycle phase shift correction in case of L2C */
        if (opt->phasshft==ISBOPT_TABLE&&isL2C(obs[i].code[1])) {
            for (j=0;j<nav->sfts.n;j++) {
                if (strcmp(nav->sfts.rectyp[j],opt->rectype[0])==0) {
                    ydif=nav->sfts.bias[j];
                }
            }
            y[1+i*nf*2]+=ydif*nav->lam[obs[i].sat-1][1];
        }
     
    }
    trace(4,"rr_=%.3f %.3f %.3f\n",rr_[0],rr_[1],rr_[2]);
    trace(4,"pos=%.9f %.9f %.3f\n",pos[0]*R2D,pos[1]*R2D,pos[2]);
    for (i=0;i<n;i++) {
        trace(4,"sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
              obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],dts[i*2],azel[i*2]*R2D,
              azel[1+i*2]*R2D);
    }
    trace(4,"y=\n"); tracemat(4,y,nf*2,n,13,3);
    
    return 1;
}
/* test valid observation data -----------------------------------------------*/
static int validobs(int *iu, int *ir, int jj, int f, int nf, double *y)
{
    int i,j;
    if (jj<0) return 0;
    i=iu[jj];
    j=ir[jj];
    if (i<0||j<0) return 0;
    /* if no phase observable, psudorange is also unusable */
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* double-differenced measurement error covariance ---------------------------*/
static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                  int nv, double *R)
{
    int i,j,k=0,b;
    
    trace(3,"ddcov   : n=%d\n",n);
    
    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {
        
        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
            R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
        }
    }
    trace(5,"R=\n"); tracemat(5,R,nv,nv,8,6);
}
/* precise tropspheric model -------------------------------------------------*/
static double prectrop(gtime_t time, const double *pos, int r,
                       const double *azel, const prcopt_t *opt, const double *x,
                       double *dtdx)
{
    double m_w=0.0,cotz,grad_n,grad_e;
    int i=IT(r,opt);
    
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
/* glonass inter-channel bias correction -------------------------------------*/
static double gloicbcorr(int sat1, int sat2, const prcopt_t *opt, double lam1,
                         double lam2, int f)
{
    double dfreq;
    
    if (f>=NFREQGLO||f>=opt->nf||!opt->exterr.ena[2]) return 0.0;
    
    dfreq=(CLIGHT/lam1-CLIGHT/lam2)/(f==0?DFRQ1_GLO:DFRQ2_GLO);
    
    return opt->exterr.gloicb[f]*0.01*dfreq; /* (m) */
}
/* test navi system (m=0:gps/qzs/sbs,1:glo,2:gal,3:bds,4:qzs) ----------------*/
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

static int check_same_fac(const obsd_t *obs, int ir[][MAXSAT])
{
    return (obs[ir[0][0]].facility==obs[ir[1][0]].facility);
}

static int get_iriu_index(int sati, int ns, const int *sat) {
    int i;
    for (i=0;i<ns;i++) {
        if(sati==sat[i]) {
            return i;
        }    
    }
    return -1;
}

static int valid_zd_sig(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, double *y,
                        int nf, int f, int j, int m, int ch, int sati, int iu[][MAXSAT], 
                        int ir[][MAXSAT], int sat[][MAXSAT], int *ff, int *ns, int *satj, 
                        int *jj, int *sysi, int *sysj, double *lami, double *lamj) {
    prcopt_t *opt=&rtk->opt;
    int code;

    *satj=obs[j].sat;
    if ((*jj=get_iriu_index(*satj,ns[ch],sat[ch]))==-1) return 0;
    if (sati==*satj) return 0;
    *sysi=rtk->ssat[sati-1].sys;
    *sysj=rtk->ssat[*satj-1].sys;
    code=rtk->ssat[*satj-1].code[f%nf];
    if (!test_sys(*sysj,m,opt->qzsmodear,code)) return 0;
    if (!validobs(iu[ch],ir[ch],*jj,f,nf,y)) return 0;

    *ff=f%nf;
    *lami=nav->lam[sati-1][*ff];
    *lamj=nav->lam[*satj-1][*ff];
    if (*lami<=0.0||*lamj<=0.0) return 0;

    return 1;
}

/* set_init_pb() is created by copying the second half of the udbias_ppp() */
static void set_init_pb(rtk_t *rtk, const nav_t *nav, const obsd_t *obs, int n, int f, int satj)
{
    int i, j, isat;
    double offset, lami, cp, pr;
    double* bias;
    bias = zeros(n, 1);

    /* estimate approximate phase-bias by phase - code */
    for (i = j = 0, offset = 0.0; i < n; i++) {
        cp = obs[i].L[f];
        pr = obs[i].P[f];
        lami = nav->lam[obs[i].sat - 1][f];
        if (cp == 0.0 || pr == 0.0 || lami == 0.0) continue;

        bias[i] = cp * lami - pr;

        /* offset = sum of (bias - phase-bias) for all valid sats in meters */
        if (rtk->x[IB(obs[i].sat, f, &rtk->opt)] != 0.0) {
            lami = nav->lam[obs[i].sat - 1][f];
            offset += bias[i] - rtk->x[IB(obs[i].sat, f, &rtk->opt)] * lami;
            j++;
        }
    }
    rtk->com_bias = j > 0 ? offset / j : 0; /* save offset for initialization below */
    /* set initial states of phase-bias */
    for (i = 0; i < n; i++) {
        isat = obs[i].sat;
        if (isat != satj) continue;
        if (bias[i] == 0.0) continue;
        trace(2, "set initial states of phase-bias(set_init_pb), sat=%2d f=%2d\n", obs[i].sat, f + 1);
        lami = nav->lam[isat - 1][f];
        initx(rtk, (bias[i] - rtk->com_bias) / lami, SQR(rtk->opt.std[0]), IB(obs[i].sat, f, &rtk->opt));
        rtk->ssat[isat - 1].lock[f] = -rtk->opt.minlock;
        rtk->ssat[isat - 1].pbreset[f] = 1;
    }
    free(bias);
}
/* double-differenced phase/code residuals -----------------------------------*/
static int ddres(rtk_t *rtk, const nav_t *nav, double *dt, const double *x,
                 const double *P, const obsd_t *obs, int n, int sat[][MAXSAT], double *y, double *e,
                 double *azel, int iu[][MAXSAT], int ir[][MAXSAT], int *ns, double *v,
                 double *H, double *R, int *vflg, int niter, int order)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im[SSR_CH_NUM];
    double *tropr[SSR_CH_NUM],*tropu[SSR_CH_NUM],*dtdxr,*dtdxu,*Ri,*Rj,lami,lamj,fi,fj,df,*Hi=NULL;
    int h,i,j,k,m,f,ff,nv=0,nb[NFREQ*6*2+1]={0},b=0,sysi,sysj,nf=NF(opt);
    static int refsat[NFREQ*2*MAXREF*6+1]={0};
    static int refsat2[NFREQ*2*6+1]={0};
    int flg=0,code,ch,sch,pch=0,samefac=0,refchgflg,ns_,sati,satj,refsatch,rch;
    int v0,v1,j0,j1,jj,vidx;
    double azeli=0.0;
    static int use_ch[NFREQ][MAXSAT], pre_use_ch[NFREQ][MAXSAT];
    double tow = time2gpst(obs[0].time, NULL);
    
    ns_=opt->l6mrg?ns[0]+ns[1]:ns[0];
    trace(2,"ddres   : tow=%.1f dt[0]=%.1f dt[1]=%.1f nx=%d ns=%d\n",tow,dt[0],dt[1],rtk->nx,ns_);
    
    bl=baseline(x,rtk->rb,dr);
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);
    memcpy(pre_use_ch, use_ch, sizeof(use_ch));
    for(f=0;f<NFREQ;f++) for(i=0;i<MAXSAT;i++) use_ch[f][i]=-1;
    
    Ri=mat(ns_*nf*2+2,1); Rj=mat(ns_*nf*2+2,1);
    dtdxu=mat(ns_,3); dtdxr=mat(ns_,3);

    for (ch=0;ch<(opt->l6mrg?SSR_CH_NUM:1);ch++) {
        tropu[ch]=mat(ns[ch],1); tropr[ch]=mat(ns[ch],1); im[ch]=mat(ns[ch],1);
    }
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].resp[j]=rtk->ssat[i].resc[j]=0.0;
    }
    /* compute factors of ionospheric and tropospheric delay */
    for (ch=0;ch<(opt->l6mrg?SSR_CH_NUM:1);ch++) {
        for (i=0;i<ns[ch];i++) {
            if (opt->ionoopt>=IONOOPT_EST) {
                im[ch][i]=(ionmapf(posu,azel+iu[ch][i]*2)+ionmapf(posr,azel+ir[ch][i]*2))/2.0;
            }
            if (opt->tropopt>=TROPOPT_EST) {
                tropu[ch][i]=prectrop(rtk->sol.time,posu,0,azel+iu[ch][i]*2,opt,x,dtdxu+i*3);
                tropr[ch][i]=prectrop(rtk->sol.time,posr,1,azel+ir[ch][i]*2,opt,x,dtdxr+i*3);
            }
        }
    }

    /* check if the same facility is used */
    if (opt->l6mrg) samefac=check_same_fac(obs,ir);

    for (m=0;m<6;m++) {/* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds,4:qzs */
    
    for (f=opt->mode>PMODE_DGPS?0:nf;f<nf*2;f++) {
        /* search reference satellite with highest elevation */
        for (h=0,sch=-1,refchgflg=0,refsatch=-1; h<((opt->l6mrg&&!samefac)?SSR_CH_NUM:1); h++) {
            azeli=0.0;j0=-1;j1=-1;
            for (i=-1,j=0,flg=0;j<n;j++) {
                sati=obs[j].sat;
                sysi=rtk->ssat[sati-1].sys;
                code=rtk->ssat[obs[j].sat-1].code[f%nf];
                if (!test_sys(sysi,m,opt->qzsmodear,code)) continue;
                if (niter>0) {
                    for (k=0;k<niter;k++) {
                        if (refsat[NFREQ*2*MAXREF*m+NFREQ*2*k+f]==sati) flg=1;
                    }
                    if (flg==1) continue;
                }
                j0=get_iriu_index(sati,ns[0],sat[0]);
                j1=get_iriu_index(sati,ns[1],sat[1]); /* only valid for 2ch */
                if (j0==-1 && j1==-1) continue;

                v0=0;v1=0;
                if (samefac || h==1) {
                    /* refsat is selected from one side ch */
                    v0=validobs(iu[0],ir[0],j0,f,nf,y);
                    v1=validobs(iu[1],ir[1],j1,f,nf,y);
                    if (!v0 && !v1) continue;
                    trace(5, "h:%d sati:%d v0:%d v1:%d ir[0][j0]:%d\n", h,sati,v0,v1,ir[0][j0==-1?j1:j0]);
                    ch=(!v0)?1:0;
                } else { /* h=0 || samefac==0 */
                    /* select common refsat between 2ch */
                    if (!(v0=validobs(iu[0],ir[0],j0,f,nf,y))) continue;
                    if (opt->l6mrg && !(v1=validobs(iu[1],ir[1],j1,f,nf,y))) continue;
                    ch=(!v0)?1:0;
                }
                if (opt->qzsmodear==2&&opt->posopt[8]&&sysi==SYS_QZS) continue;

                vidx=v0?j0:j1;
                if (i<0||azel[1+iu[ch][vidx]*2]>=azeli){
                    i=j;
                    azeli= azel[1+iu[ch][vidx]*2];
                    if (samefac) refsatch=ch;
                    if (h==1) sch=ch;
                    trace(5, "refsatroop m:%d i:%d sat:%d azel:%f\n", m, i, sati, azeli);
                }
            }
            if (i>=0) break;
        }
        if (i<0) continue;
        sati = obs[i].sat;
        if (f<nf) use_ch[f][sati]=2; /* use_ch=2 means ref sat */
        if (f<nf && refsat2[NFREQ*2*m+NFREQ*2+f]!=sati && refsat2[NFREQ*2*m+NFREQ*2+f]!=0) {
            trace(2,"refsat changed2 %s L%d sat=%2d -> %2d\n",time_str(rtk->sol.time,0),
                    f%nf+1,refsat2[NFREQ*2*m+NFREQ*2+f], sati);
            refchgflg =1;
        }
        refsat[NFREQ*2*MAXREF*m+NFREQ*2*niter+f]=sati;
        refsat2[NFREQ*2*m+NFREQ*2+f]=sati;
        if (niter>0) {
           trace(2,"refsat changed %s L%d sat=%2d -> %2d\n",time_str(rtk->sol.time,0),
                   f%nf+1,refsat[NFREQ*MAXREF*2*m+NFREQ*2*(niter-1)+f],
                  refsat[NFREQ*MAXREF*2*m+NFREQ*2*niter+f]);
        }
        if (f<nf) trace(4, "refsat m:%d L%d %s: %d niter:%d refsatch:%d\n",m,f+1,time_str(rtk->sol.time,0),sati,niter,refsatch);

        /* make double difference */
        if (opt->l6mrg) pch=(ns[0]>=ns[1])?0:1;
        for (j=0;j<n;j++) {
            if (samefac) { rch=refsatch; } else { rch=(sch!=-1)?sch:pch; }
            i=get_iriu_index(sati,ns[rch],sat[rch]);
            if (!opt->l6mrg) {
                ch=0; 
                if (!valid_zd_sig(rtk,obs,nav,y,nf,f,j,m,ch,sati,iu,ir,sat,&ff,ns,&satj,&jj,&sysi,&sysj,&lami,&lamj)) continue;
            } else if (sch==-1) {
                /* priority channel selection */
                ch=pch;
                if (!valid_zd_sig(rtk,obs,nav,y,nf,f,j,m,ch,sati,iu,ir,sat,&ff,ns,&satj,&jj,&sysi,&sysj,&lami,&lamj)) {
                    ch=(ch==1)?0:1;
                    if (!samefac) { i=get_iriu_index(sati,ns[ch],sat[ch]); rch=ch; }
                    if (!valid_zd_sig(rtk,obs,nav,y,nf,f,j,m,ch,sati,iu,ir,sat,&ff,ns,&satj,&jj,&sysi,&sysj,&lami,&lamj)) continue;
                }
            } else {
                trace(3,"even though it is a clas 2ch case, refsat is selected from one side ch.(ch=%d).\n", sch);
                ch=sch;
                if (!valid_zd_sig(rtk,obs,nav,y,nf,f,j,m,ch,sati,iu,ir,sat,&ff,ns,&satj,&jj,&sysi,&sysj,&lami,&lamj)) continue;
            }
            if (f<nf) use_ch[f][satj]=ch;
            if (order==FST && f<nf && opt->l6mrg && !samefac) {
                if ((pre_use_ch[f][sati]!=ch && refchgflg) || (pre_use_ch[f][satj]!=ch)) {
                    set_init_pb(rtk, nav, obs, n, f, satj);
                }
            }

            if (H) Hi=H+nv*rtk->nx;
            
            /* double-differenced residual */
            v[nv]=(y[f+iu[rch][i]*nf*2]-y[f+ir[rch][i]*nf*2])-
                  (y[f+iu[ch][jj]*nf*2]-y[f+ir[ch][jj]*nf*2]);

            if (order==FST) {
                trace(5, "dd_res sati:%d satj:%d rch:%d ch:%d v:%f\n", sati, satj, rch, ch, v[nv]);
            }
    
            /* partial derivatives by rover position */
            if (H) {
                for (k=0;k<3;k++) {
                    Hi[k]=-e[k+iu[rch][i]*3]+e[k+iu[ch][jj]*3];
                }
            }
            /* double-differenced ionospheric delay term */
            if (opt->ionoopt==IONOOPT_EST||opt->ionoopt==IONOOPT_EST_ADPT) {
                fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
                didxi=(f<nf?-1.0:1.0)*fi*fi*im[rch][i];
                didxj=(f<nf?-1.0:1.0)*fj*fj*im[ch][jj];
                v[nv]-=didxi*x[II(sati,opt)]-didxj*x[II(satj,opt)];
                if (H) {
                    Hi[II(sati,opt)]= didxi;
                    Hi[II(satj,opt)]=-didxj;
                }
            }
            /* double-differenced tropospheric delay term */
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                v[nv]-=(tropu[rch][i]-tropu[ch][jj])-(tropr[rch][i]-tropr[ch][jj]);
                for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {
                    if (!H) continue;
                    Hi[IT(0,opt)+k]= (dtdxu[k+i*3]-dtdxu[k+jj*3]);
                    Hi[IT(1,opt)+k]=-(dtdxr[k+i*3]-dtdxr[k+jj*3]);
                }
            }
            /* double-differenced phase-bias term */
            if (f<nf) {
                if (opt->ionoopt!=IONOOPT_IFLC) {
                    v[nv]-=lami*x[IB(sati,f,opt)]-lamj*x[IB(satj,f,opt)];
                    if (H) {
                        Hi[IB(sati,f,opt)]= lami;
                        Hi[IB(satj,f,opt)]=-lamj;
                    }
                }
                else {
                    v[nv]-=x[IB(sati,f,opt)]-x[IB(satj,f,opt)];
                    if (H) {
                        Hi[IB(sati,f,opt)]= 1.0;
                        Hi[IB(satj,f,opt)]=-1.0;
                    }
                }
            }
            /* glonass receiver h/w bias term */
            if (rtk->opt.glomodear==2&&sysi==SYS_GLO&&sysj==SYS_GLO&&ff<NFREQGLO) {
                df=(CLIGHT/lami-CLIGHT/lamj)/1E6; /* freq-difference (MHz) */
                v[nv]-=df*x[IL(ff,opt)];
                if (H) Hi[IL(ff,opt)]=df;
            }
            /* glonass interchannel bias correction */
            else if (sysi==SYS_GLO&&sysj==SYS_GLO) {
                
                v[nv]-=gloicbcorr(sati,satj,&rtk->opt,lami,lamj,f);
            }
#if 0
            /* test innovation */
            if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                if (f<nf) {
                    rtk->ssat[sat[i]-1].rejc[f]++;
                    rtk->ssat[sat[j]-1].rejc[f]++;
                }
                errmsg(rtk,"outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                       sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                continue;
            }
#endif

#if 1
            if (f<nf) rtk->ssat[satj-1].resc[f   ]=v[nv];
            else      rtk->ssat[satj-1].resp[f-nf]=v[nv];
#else
            if (f<nf) rtk->ssat[sat[j]-1].resc[f   ]=-v[nv];
            else      rtk->ssat[sat[j]-1].resp[f-nf]=-v[nv];
#endif

            /* single-differenced measurement error variances */
            Ri[nv]=varerr(sati,sysi,azel[1+iu[rch][i]*2],bl,dt[rch],f,opt);
            Rj[nv]=varerr(satj,sysj,azel[1+iu[ch][jj]*2],bl,dt[ch],f,opt);
            if (f==1) {
                Ri[nv]*=SQR(2.55/1.55);
                Rj[nv]*=SQR(2.55/1.55);
            }            
            /* set valid data flags */
            if (f<nf) rtk->ssat[sati-1].vsat[f]=rtk->ssat[satj-1].vsat[f]=1;
            
            trace(4,"sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",sati,
                  satj,f<nf?"L":"P",f%nf+1,v[nv],Ri[nv],Rj[nv]);
            
            vflg[nv++]=(sati<<16)|(satj<<8)|((f<nf?0:1)<<4)|(f%nf);
            nb[b]++;
        }
#if 0
        if(f<nf){
            if (rtk->ssat[sat[i]-1].rejc[f]<0) rtk->ssat[sat[i]-1].rejc[f]=0;
        }
#endif
#if 0 /* residuals referenced to reference satellite (2.4.2 p11) */
        /* restore single-differenced residuals assuming sum equal zero */
        if (f<nf) {
            for (j=0,s=0.0;j<MAXSAT;j++) s+=rtk->ssat[j].resc[f];
            s/=nb[b]+1;
            for (j=0;j<MAXSAT;j++) {
                if (j==sat[i]-1||rtk->ssat[j].resc[f]!=0.0) rtk->ssat[j].resc[f]-=s;
            }
        }
        else {
            for (j=0,s=0.0;j<MAXSAT;j++) s+=rtk->ssat[j].resp[f-nf];
            s/=nb[b]+1;
            for (j=0;j<MAXSAT;j++) {
                if (j==sat[i]-1||rtk->ssat[j].resp[f-nf]!=0.0)
                    rtk->ssat[j].resp[f-nf]-=s;
            }
        }
#endif
        b++;
        }
    }
    /* end of system loop */
    
    /* double-differenced measurement error covariance */
    ddcov(nb,b,Ri,Rj,nv,R);
    
    free(Ri); free(Rj); free(dtdxu); free(dtdxr);
    for (ch=0;ch<(opt->l6mrg?SSR_CH_NUM:1);ch++) {
        free(tropu[ch]); free(tropr[ch]); free(im[ch]);
    }
    
    return nv;
}
/* time-interpolation of residuals (for post-mission) ------------------------*/
static double intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                      rtk_t *rtk, int ch, double *y)
{
    static obsd_t obsb[SSR_CH_NUM][MAXOBS];
    static double yb[SSR_CH_NUM][MAXOBS*NFREQ*2],rs[SSR_CH_NUM][MAXOBS*6],dts[SSR_CH_NUM][MAXOBS*2],var[SSR_CH_NUM][MAXOBS];
    static double e[SSR_CH_NUM][MAXOBS*3],azel[SSR_CH_NUM][MAXOBS*2];
    static int nb[SSR_CH_NUM]={0},svh[SSR_CH_NUM][MAXOBS*2];
    prcopt_t *opt=&rtk->opt;
    double tt=timediff(time,obs[0].time),ttb,*p,*q;
    int i,j,k,nf=NF(opt);
    
    trace(3,"intpres : n=%d tt=%.1f\n",n,tt);
    
    if (nb[ch]==0||fabs(tt)<DTTOL) {
        nb[ch]=n; for (i=0;i<n;i++) obsb[ch][i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[ch][0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;
    
    satposs(time,obsb[ch],nb[ch],nav,opt->sateph,rs[ch],dts[ch],var[ch],svh[ch]);
    
    if (!zdres(rtk,1,obsb[ch],nb[ch],rs[ch],dts[ch],svh[ch],nav,rtk->rb,opt,1,yb[ch],e[ch],azel[ch])) {
        return tt;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nb[ch];j++) if (obsb[ch][j].sat==obs[i].sat) break;
        if (j>=nb[ch]) continue;
        for (k=0,p=y+i*nf*2,q=yb[ch]+j*nf*2;k<nf*2;k++,p++,q++) {
            if (*p==0.0||*q==0.0) *p=0.0; else *p=(ttb*(*p)-tt*(*q))/(ttb-tt);
        }
    }
    return fabs(ttb)>fabs(tt)?ttb:tt;
}
/* single to double-difference transformation matrix (D') --------------------*/
static int ddmat(rtk_t *rtk, double *D)
{
    int i,j,k,m,f,nb=0,nx=rtk->nx,na=rtk->na,nf=NF(&rtk->opt),code;
    
    trace(3,"ddmat   :\n");
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].fix[j]=0;
    }
    for (i=0;i<na;i++) D[i+i*nx]=1.0;
    
    for (m=0;m<6;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds,4:qzs */
        
        if (m==1&&rtk->opt.glomodear==0) continue;
        if (m==3&&rtk->opt.bdsmodear==0) continue;
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
    
    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i];
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i];
    
    for (m=0;m<6;m++) for (f=0;f<nf;f++) {
        
        for (n=i=0;i<MAXSAT;i++) {
            code=rtk->ssat[i].code[f];
            if (!test_sys(rtk->ssat[i].sys,m,rtk->opt.qzsmodear,code)||
                          rtk->ssat[i].fix[f]!=2) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
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
    int i,n,m,f,info,index[MAXSAT],nb=rtk->nx-rtk->na,nv=0,nf=NF(&rtk->opt);
    int code;

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
            errmsg(rtk,"filter error (info=%d)\n",info);
        }
        free(R);
    }
    free(v); free(H);
}
/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,ny,nb,info,nx=rtk->nx,na=rtk->na;
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*QQ,s[2];
    double *Ds,*Pac,*Pcc,*D_P;
    int nc;

    trace(3,"resamb_LAMBDA : nx=%d\n",nx);
    
    rtk->sol.ratio=0.0;
    
    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D=zeros(nx,nx);
    if ((nb=ddmat(rtk,D))<=0) {
        errmsg(rtk,"no valid double-difference\n");
        free(D);
        return 0;
    }
    ny=na+nb; y=mat(ny,1); Qy=mat(ny,ny); DP=mat(ny,nx);
    b=mat(nb,2); db=mat(nb,1); Qb=mat(nb,nb); Qab=mat(na,nb); QQ=mat(na,nb);
    
    nc=nx-na;Ds=mat(nc,nb); Pac=mat(na,nc); Pcc=mat(nc,nc);D_P=mat(nb,nc);
    for (i=0;i<nc;i++) for (j=0;j<nb;j++) Ds[i+j*nc]=D[na+i+(na+j)*nx];
    for (i=0;i<na;i++) for (j=0;j<nc;j++) Pac[i+j*na]=rtk->P[   i+(na+j)*nx];
    for (i=0;i<nc;i++) for (j=0;j<nc;j++) Pcc[i+j*nc]=rtk->P[na+i+(na+j)*nx];
    matmul("TN",ny, 1,nx,1.0,D ,rtk->x,0.0,y );
    matmul("NN",na,nb,nc,1.0,Pac,Ds ,0.0,Qab);
    matmul("TN",nb,nc,nc,1.0,Ds ,Pcc,0.0,D_P);
    matmul("NN",nb,nb,nc,1.0,D_P,Ds,0.0,Qb);
    free(Ds);free(Pac);free(Pcc);free(D_P);

    trace(4,"N(0)="); tracemat(4,y+na,1,nb,10,3);
    
    /* lambda/mlambda integer least-square estimation */
    if (!(info=lambda(nb,2,y+na,Qb,b,s))) {
        
        trace(4,"N(1)="); tracemat(4,b   ,1,nb,10,3);
        trace(4,"N(2)="); tracemat(4,b+nb,1,nb,10,3);
        
        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;
        
        /* validation by popular ratio-test */
        if (s[0]<=0.0||s[1]/s[0]>=qfrtk[opt->alphaar][nb-1]) {
            
            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                rtk->xa[i]=rtk->x[i];
                for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i]=b[i];
                y[na+i]-=b[i];
            }
            if (!matinv(Qb,nb)) {
                matmul("NN",nb,1,nb, 1.0,Qb ,y+na,0.0,db);
                matmul("NN",na,1,nb,-1.0,Qab,db  ,1.0,rtk->xa);
                
                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN",na,nb,nb, 1.0,Qab,Qb ,0.0,QQ);
                matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,rtk->Pa);
                
                trace(3,"resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                      nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);
                
                /* restore single-differenced ambiguity */
                restamb(rtk,bias,nb,xa);
            }
            else nb=0;
        }
        else { /* validation failed */
#if 0
            errmsg(rtk,"ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                   nb,s[1]/s[0],s[0],s[1]);
#endif
            nb=0;
        }
    }
    else {
        errmsg(rtk,"lambda error (info=%d)\n",info);
    }
    free(D); free(y); free(Qy); free(DP);
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
static void resetfilterflag(nav_t *nav)
{
    nav->filreset = FALSE;
    clearsatcorr();
}

static void check_vrs_facility(nav_t *nav, const obsd_t *obs, int nu, int *nr, int l6mrg)
{
    static int savefacility[SSR_CH_NUM] = {-1, -1};
    int ch,idx;

    for (ch=0;ch<(l6mrg?SSR_CH_NUM:1);ch++) {
        idx=(ch==0)?nu:nu+nr[0];
        if (savefacility[ch] != -1 && savefacility[ch] != obs[idx].facility) {
            trace(2, "VRS facility changed[ch:%d], %d ---> %d\n", ch, savefacility[ch], obs[idx].facility);
            nav->filreset = TRUE;
        }
        nav->facility[ch] = obs[idx].facility;
        savefacility[ch]  = obs[idx].facility;
    }
}

/* relative positioning ------------------------------------------------------*/
extern int relposvrs(rtk_t *rtk, const obsd_t *obs, int nu, int *nr, nav_t *nav)
{
    static gtime_t regularly = {-1, 0.0};
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*v,*H,*R,*xp,*Pp,*Qp,*xa,*bias,dt[SSR_CH_NUM];
    int i,j,k=0,f,l,n,ns[SSR_CH_NUM]={0},ny,nv=0,sat[SSR_CH_NUM][MAXSAT],iu[SSR_CH_NUM][MAXSAT],ir[SSR_CH_NUM][MAXSAT],sati;
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;
    int nf=opt->ionoopt==IONOOPT_IFLC?1:opt->nf;
    int ch, ofst;
    double pdop;
    int nb,isat[MAXSAT];
    double el[MAXSAT],dist;
    static int resetcnt=0,cntdiffp=0,np;
    double tow = time2gpst(obs[0].time, NULL);
    
    trace(2,"relposvrs  : time=%s nx=%d nu=%d nr[0]=%d nr[1]=%d\n",time_str(obs[0].time,0),rtk->nx,nu,nr[0],nr[1]);

    n=opt->l6mrg?nu+nr[0]+nr[1]:nu+nr[0];

    check_vrs_facility(nav,obs,nu,nr,opt->l6mrg);
    
    regularly = (regularly.time == -1 ? obs[0].time: regularly);
    if (opt->regularly != 0 && timediff(obs[0].time, regularly) >= (double)opt->regularly) {
        trace(1, "relposvrs(): regularly reset filter, tow=%.1f\n", tow);
        for (i=0;i<rtk->nx;i++) rtk->x[i]=0.0;
        rtk->sol.stat=SOLQ_SINGLE;
        regularly = obs[0].time;
        resetfilterflag(nav);
        resetcnt=0;
        return 0;
    }
    
    if (nav->filreset == TRUE) {
        trace(1, "relposvrs(): reset filter, tow=%.1f\n", time2gpst(obs[0].time, NULL));
        for (i=NP(opt);i<rtk->nx;i++) rtk->x[i]=0.0;
        rtk->sol.stat=SOLQ_FLOAT;
        resetfilterflag(nav);
        stat = SOLQ_FLOAT;
        resetcnt=0;
    }
    ++resetcnt;

    /* temporal update of states */
    udstate(rtk,obs,nu,nav);

    dt[0]=timediff(time,obs[nu].time);
    if (opt->l6mrg) dt[1]=timediff(time,obs[nu+nr[0]].time);
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); y=mat(nf*2,n); e=mat(3,n);
    azel=zeros(2,n);
    
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat[j]=rtk->ssat[i].snr[j]=0;
    }
    /* satellite positions/clocks */
    satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);
    
    for (ch=0;ch<(opt->l6mrg?SSR_CH_NUM:1);ch++) {
        ofst=ch==1?nu+nr[0]:nu;
        /* undifferenced residuals for base station */
        trace(5,"zdres: tow:%f ch:%d\n", time2gpst(obs[0].time, NULL), ch);
        if (!zdres(rtk,1,obs+ofst,nr[ch],rs+ofst*6,dts+ofst*2,
                    svh+ofst,nav,rtk->rb,opt,1,y+ofst*nf*2,
                    e+ofst*3,azel+ofst*2)) {
            errmsg(rtk,"initial base station position error\n");
        
            free(rs); free(dts); free(var); free(y); free(e); free(azel);
            return 0;
        }
        /* time-interpolation of residuals (for post-processing) */
        if (opt->intpref) {
            dt[ch]=intpres(time,obs+ofst,nr[ch],nav,rtk,ch,y+ofst*nf*2);
        }
        /* select common satellites between rover and base-station */
        ns[ch]=selsat(obs,azel,nu,nr,opt,sat[ch],iu[ch],ir[ch],ch);
    }
    if (ns[0]+ns[1]==0) {
        errmsg(rtk,"no common satellite\n");
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    
    trace(4,"x(0)="); tracemat(4,rtk->x,1,NR(opt),13,4);
    
    xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx); xa=mat(rtk->nx,1);
    matcpy(xp,rtk->x,rtk->nx,1);
    
    ny=(ns[0]+ns[1])*nf*2+2;
    v=mat(ny,1); H=zeros(rtk->nx,ny); R=mat(ny,ny); bias=mat(rtk->nx,1);
    Qp=zeros(rtk->nx,rtk->nx);
    
    ch=0;
    for (i=0;i<rtk->opt.niter;i++) {
        for (k=0;k<MAXREF;k++) {
            nv=zdres(rtk,0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel);
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
            if ((nv=ddres(rtk,nav,dt,xp,Pp,obs,nu,sat,y,e,azel,iu,ir,ns,v,H,R,vflg,k,FST))<=0) {
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
            if (zdres(rtk,0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel)) {
                nv=ddres(rtk,nav,dt,xp,Pp,obs,nu,sat,y,e,azel,iu,ir,ns,v,H,R,vflg,k,SND);
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

    /* update process noise of ionoshepre delay */
    if (rtk->opt.ionoopt==IONOOPT_EST_ADPT) {
        for (i=0;i<nu;i++) {
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
        for (i=0;i<nu;i++) {
            sati = obs[i].sat;
            for (f=0;f<nf;f++) {
                if (!rtk->ssat[sati-1].vsat[f]) continue;
                rtk->ssat[sati-1].lock[f]++;
                rtk->ssat[sati-1].outc[f]=0;
                if (f==0) rtk->sol.ns++; /* valid satellite count by L1 */
            }
        }
        /* lack of valid satellites */
        if (rtk->sol.ns<4) stat=SOLQ_NONE;
    }
    
    /* keep symmetry of covariance matrix */
    for (i=0;i<rtk->nx;i++) for (j=0;j<i;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=(rtk->P[i+j*rtk->nx]+rtk->P[j+i*rtk->nx])*0.5;
    }

    /* resolve integer ambiguity by LAMBDA */
    pdop=calpdop(rtk,H,azel,vflg,obs,rtk->opt.elmaskar,nu,nv);
    if (stat!=SOLQ_NONE&&(rtk->opt.maxpdopar==0.0||pdop<=rtk->opt.maxpdopar)) {
        if ((nb=resamb_LAMBDA(rtk,bias,xa))<=1&&rtk->opt.posopt[6]) {
            int exsat=0;
            double ratio=0.0;
            elsort(rtk,obs,isat,el,nu);
            for (l=0;l<rtk->opt.armaxdelsat;l++) {
                for (i=0;i<nu;i++) {
                    sati = isat[i];
                    if (el[i]<=rtk->opt.elmaskar) continue;
                    for (f=j=0;f<nf;f++) {
                        if (!rtk->ssat[sati-1].vsat[f]||rtk->ssat[sati-1].lock[f]<=0) j++;
                    }
                    if (j>1) continue;
                    for (f=0;f<nf;f++) rtk->ssat[sati-1].vsat[f]=0;
                    pdop=calpdop(rtk,H,azel,vflg,obs,rtk->opt.elmaskar,nu,nv);
                    if ((rtk->opt.maxpdopar==0.0||pdop<=rtk->opt.maxpdopar)&&
                                        (nb=resamb_LAMBDA(rtk,bias,xa))>1) {
                        trace(2,"PAR OK l=%d i=%2d nb=%2d excluded sat=%2d el=%4.1f ratio=%5.2f\n",
                              l,i,nb,sati,el[i]*R2D,rtk->sol.ratio);
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
        
        if (nb>1&&zdres(rtk,0,obs,nu,rs,dts,svh,nav,rtk->opt.mode==PMODE_FIXED?rtk->opt.ru:xa,opt,0,y,e,azel)) {
            /* post-fix residuals for fixed solution */
            nv=ddres(rtk,nav,dt,xa,NULL,obs,nu,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg,k,TRD);
            /* validation of fixed solution */
            filter2(rtk,xp,Pp,NULL,H,v,R,rtk->nx,nv,vflg,2);
            
            if (rtk->sol.chisq<opt->maxinno[3]) {
                /* hold integer ambiguity */
                pdop=calpdop(rtk,H,azel,vflg,obs,rtk->opt.elmaskhold,nu,nv);
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
            rtk->sol.rr[i]=rtk->xa[i];
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[1+2*rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];
        rtk->nfix=0;
    }
    for (i=0;i<n;i++) for (j=0;j<nf;j++) {
        if (obs[i].L[j]==0.0) continue;
        rtk->ssat[obs[i].sat-1].pt[obs[i].rcv-1][j]=obs[i].time;
        rtk->ssat[obs[i].sat-1].ph[obs[i].rcv-1][j]=obs[i].L[j];
    }
    for (i=0;i<ns[ch];i++) for (j=0;j<nf;j++) {
        
        /* output snr of rover receiver */
        rtk->ssat[sat[ch][i]-1].snr[j]=obs[iu[ch][i]].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<nf;j++) {
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
        if (rtk->ssat[i].slip[j]&1) rtk->ssat[i].slipc[j]++;
    }

    free(rs); free(dts); free(var); free(y); free(e); free(azel);
    free(xp); free(Pp);  free(xa);  free(v); free(H); free(R); free(bias);
    free(Qp);
    
    if (stat!=SOLQ_NONE) rtk->sol.stat=stat;

    if (stat==SOLQ_NONE) {
        for (i=0;i<rtk->nx;i++) rtk->x[i]=0.0;
    }
    
    return stat!=SOLQ_NONE;
}

