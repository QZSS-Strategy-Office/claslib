/*------------------------------------------------------------------------------
* cssr2osr.c : cssr to osr converter
*
*          Copyright (C) 2015- by Mitsubishi Electric Corporation, All rights reserved.
*          Copyright (C) 2010- by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL use IERS tide model
*
*          
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
* history : 2018/03/29 1.0  new as claslib
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "cssr2osr.h"

#ifdef ENA_PPP_RTK
extern FILE *fp_osr;
#endif
/* global variables ----------------------------------------------------------*/
#ifdef CSSR2OSR_VRS
static unsigned char obsfreqs[]={ /* 1:L1,2:L2,3:L5,4:L6,5:L7,6:L8,7:L3 */
    
    0, 1, 1, 1, 1,  1, 1, 1, 1, 1, /*  0- 9 */
    1, 1, 1, 1, 2,  2, 2, 2, 2, 2, /* 10-19 */
    2, 2, 2, 2, 3,  3, 3, 5, 5, 5, /* 20-29 */
    4, 4, 4, 4, 4,  4, 4, 6, 6, 6, /* 30-39 */
    2, 2, 4, 4, 3,  3, 3, 1, 1, 0  /* 40-49 */
};
#endif

/* precise tropospheric model ------------------------------------------------*/
extern double prectrop(gtime_t time, const double *pos, const double *azel,
                       const prcopt_t *opt, const double zwd, double ztd)
{
    double m_d,m_w,tdvd,twvd,gh;

    gh=geoidh(pos);
    if (get_stTv(time,pos[0],pos[2],gh,&tdvd,&twvd)) return 0.0;

    /* mapping function */
    m_d=tropmapf(time,pos,azel,&m_w);

    return m_d*tdvd*ztd+m_w*twvd*zwd;

}

/* phase and code residuals --------------------------------------------------*/
#ifndef CSSR2OSR_VRS
extern int zdres(const obsd_t *obs,
#else
extern int zdres(obsd_t *obs,
#endif
                 int n, const double *rs, const double *dts,
                 const double *vare, const int *svh, nav_t *nav,
                 double *x, double *y, 
                 double *e, double *azel, rtk_t *rtk,
                 int osrlog, double *cpc, gtime_t *pt0, grid_t *grid,
                 ssat_t *ssat, prcopt_t *opt, sol_t *sol, osrd_t *osr)
{
    double r,rr[3],disp[3],pos[3],meas[NFREQ*2],modl[NFREQ*2],zwd,ztd,quality;
    double *lam,fi;
    int i,j,f,qj,sat,sys,nv=0,brk,nf=opt->nf;
    int tbrk=0;
    double isb[NFREQ][2]={0},ydif=0.0;
    int nsig,k=0;
    double isb_by_prn[NFREQ];
    double tow;
    static double pbias_ofst[MAXSAT*(NFREQ+NEXOBS)]={0};
    double dcpc;
#ifdef ENA_PPP_RTK
    osrd_t osrtmp[MAXOBS]={{{0}}};
    osr = osrtmp;
#else
    static double cpctmp[MAXSAT*(NFREQ+NEXOBS)]={0};
    static gtime_t pt0tmp[MAXSAT*(NFREQ+NEXOBS)]={0};
    pt0=pt0tmp;
    cpc=cpctmp;
#endif

#ifndef CSSR2OSR_VRS    
    if (n<=0) return 0;
#endif

#ifdef ENA_PPP_RTK
    for (i=0;i<n*nf*2;i++) y[i] = 0.0;
#endif

    trace(3,"zdres : n=%d\n",n);

    for (i=0;i<3;i++) rr[i]=x[i];

    if (norm(rr,3)<=0.0) return 0;

    sol->network = grid->network;
    ecef2pos(rr,pos);
    if (!trop_grid_data(nav, grid->index, obs[0].time, grid->num, grid->weight, grid->Gmat, grid->Emat, &zwd, &ztd, &quality,&tbrk)) { /* zenith trop delay */
        trace(2,"trop correction error: time=%s \n", time_str(obs[0].time,2));
        return 0;
    }

    /* earth tides correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rr,opt->tidecorr,&nav->erp,opt->odisp[0],
                 nav->oload[grid->network-1],grid->network,grid->num,grid->index,grid->Gmat,grid->Emat,grid->weight,disp);
        for (i=0;i<3;i++) rr[i]+=disp[i];
        ecef2pos(rr,pos);
    }

#ifdef ENA_PPP_RTK
    /* check invalid grid */
    nav->invtrop=tbrk;
#endif

    for (i=0,k=0;i<n&&i<MAXOBS;i++) {
        int iodeflag = FALSE, prn;
        double tmp_r, tmp_dts;
        sat=obs[i].sat;
        lam=nav->lam[sat-1];
        sys=satsys(sat,&prn);
        osr[i].sat=sat;
        
        for (j=0;j<nf*2;j++) meas[j]=modl[j]=0.0;

        /* geometric distance/azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr,e+i*3))<=0.0) {
            continue;
        }
        if (satazel(pos,e+i*3,azel+i*2)<opt->elmin) {
            trace(3,"elevation mask rejection: sat=%2d,el=%4.1f\n",sat,azel[i*2+1]*R2D);
            continue;
        }

        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;

        /* shapiro time delay correction */
        if (opt->posopt[7]) osr[i].relatv=shapiro(rs+i*6,rr);

        /* tropospheric delay correction */
        osr[i].trop=prectrop(obs[i].time,pos,azel+i*2,opt,zwd,ztd);

        /* phase windup correction */
        windupcorr(sol->time,rs+i*6,rr,&ssat[sat-1].phw);

        /* ionosphere , satellite code/phase bias and phase windup corrected measurements */
        if (!corrmeas(obs+i, nav, pos, azel+i*2, opt,
                      grid->index,grid->num,grid->weight,grid->Gmat,grid->Emat,ssat[sat-1],&brk,osr+i,ssat[sat-1].pbreset,rtk->sisadjust)) {
            continue;
        }

#if 0
        /* low reference ambiguity fix ?  */
        if (opt->threexclssr>0.0&&(nav->ssr[sat-1].fixratio<opt->threexclssr)) {
            trace(2,"ssr rejection due to low reference ambiguity fix ratio:sys=%2d sat=%2d thres=%3.2f ratio=%3.2f\n",
                sys,sat,opt->threexclssr,nav->ssr[sat-1].fixratio);
            continue;
        }
#endif

#ifndef CSSR2OSR_VRS
        for (j=0;j<nf;j++) {
            f=j;
#else
        nsig=nav->ssr[sat-1].nsig;
        for (j=0;j<nsig;j++) {
            f=obsfreqs[nav->ssr[sat-1].smode[j]];f=f>0?f-1:0;
#endif
            fi=lam[f]/lam[0];
            if (osr[i].pbias[j]==CSSRINVALID||osr[i].cbias[j]==CSSRINVALID) {
#ifndef ENA_PPP_RTK
                osr[i].cbias[j]=osr[i].cbias[j]==CSSRINVALID?0.0:osr[i].cbias[j];
                osr[i].pbias[j]=osr[i].pbias[j]==CSSRINVALID?0.0:osr[i].pbias[j];
#endif
                continue;
            }
            osr[i].PRC[j]=osr[i].trop+osr[i].relatv+osr[i].antr[f]+fi*fi*FREQ2/FREQ1*osr[i].iono
                            +osr[i].cbias[j];
            osr[i].CPC[j]=osr[i].trop+osr[i].relatv+osr[i].antr[f]-fi*fi*FREQ2/FREQ1*osr[i].iono
                            +osr[i].pbias[j]+osr[i].wupL[f]+osr[i].compL[j];
            if (opt->posopt[11] != 0) {
                if (!opt->posopt[9]) {
                    osr[i].CPC[j] = adjust_cpc(obs[i].time, sat, &nav->ssr[sat-1], f, osr[i].CPC[j], &osr[i].sis, &iodeflag);
                    osr[i].PRC[j] = adjust_prc(obs[i].time, sat, &nav->ssr[sat-1], f, osr[i].PRC[j], &osr[i].sis, &iodeflag);
                }
            } else {
                if (rtk->sisadjust) {
                    osr[i].CPC[j] = adjust_cpc(obs[i].time, sat, &nav->ssr[sat-1], f, osr[i].CPC[j], &osr[i].sis, &iodeflag);
                    osr[i].PRC[j] = adjust_prc(obs[i].time, sat, &nav->ssr[sat-1], f, osr[i].PRC[j], &osr[i].sis, &iodeflag);
                }
            }
            if (iodeflag == TRUE) {
                if (j == 0) {
                    /* transmission time by satellite clock */
                    gtime_t dtsat = timeadd(obs[i].time,-obs[i].P[0]/CLIGHT);
                    tmp_r = r;
                    adjust_r_dts(&tmp_r, &tmp_dts, obs[i].time, sat, nav, rr, dtsat);
                    trace(4, "adjust_r_dts: tow=%.1f, sat=%d, r=%.3f --> %.3f, dts=%.16f --> %.16f\n",
                          time2gpst(obs[i].time, NULL), sat, r, tmp_r, dts[i*2], tmp_dts);
                }

                modl[j]   =tmp_r-CLIGHT*tmp_dts+osr[i].CPC[j];
                modl[nf+j]=tmp_r-CLIGHT*tmp_dts+osr[i].PRC[j];
                osr[i].CPC[j]+=tmp_r-CLIGHT*tmp_dts-(r-CLIGHT*dts[i*2]);
                osr[i].PRC[j]+=tmp_r-CLIGHT*tmp_dts-(r-CLIGHT*dts[i*2]);
            } else {
                modl[j]   =r-CLIGHT*dts[i*2]+osr[i].CPC[j];
                modl[nf+j]=r-CLIGHT*dts[i*2]+osr[i].PRC[j];
            }
            osr[i].p[j] = modl[nf+j];
            osr[i].c[j] = modl[j];

#ifdef CSSR2OSR_VRS
            obs[k].P[j] = osr[i].p[j];
            obs[k].L[j] = osr[i].c[j]/lam[f];
            obs[k].LLI[j] = nav->filreset == TRUE?1:0;
#endif

            /* repair cycle slip of pbias */
            tow = time2gpst(obs[i].time, NULL);
            getorbitclock(tow, sat, &osr[i].orb, &osr[i].clk);
            if (nav->filreset == FALSE&&
                fabs(timediff(nav->ssr[sat-1].t0[5],pt0[sat-1]))>0.0&&
                fabs(timediff(nav->ssr[sat-1].t0[5],pt0[sat-1]))<120.0) {
                dcpc=osr[i].orb-osr[i].clk+osr[i].CPC[j]-cpc[j*MAXSAT+sat-1];
                if (dcpc>=95.0*lam[f]&&dcpc<105.0*lam[f]) {
                    pbias_ofst[j*MAXSAT+sat-1]-=100.0;

#ifdef ENA_PPP_RTK
                    x[IB(sat,f,opt)]-=100.0;
#endif
                    trace(2,"pbias slip detected t=%s sat=%2d f=%1d dcpc[cycle]=%.1f\n",time_str(obs[i].time,0),sat,f,dcpc/lam[f]);
                } else if (dcpc<=-95.0*lam[f]&&dcpc>-105.0*lam[f]) {
                    pbias_ofst[j*MAXSAT+sat-1]+=100.0;
#ifdef ENA_PPP_RTK
                    x[IB(sat,f,opt)]+=100.0;
#endif
                    trace(2,"pbias slip detected t=%s sat=%2d f=%1d dcpc[cycle]=%.1f\n",time_str(obs[i].time,0),sat,f,dcpc/lam[f]);
                }
            }else{
                pbias_ofst[j*MAXSAT+sat-1]=0.0;
            }
#ifdef CSSR2OSR_VRS
            obs[k].L[j] +=pbias_ofst[j*MAXSAT+sat-1];
#endif
            cpc[j*MAXSAT+sat-1]=osr[i].orb-osr[i].clk+osr[i].CPC[j];
        }
        /* isb correction */
        if (opt->isb==ISBOPT_TABLE) {
            chk_isb(satsys(sat,NULL), opt, &nav->stas[obs->rcv-1], isb);
            for (f=0;f<nf;f++) meas[f]=meas[f]-isb[f][0]; /* L */

            if (opt->isbbyprn && chk_isb_by_prn(sat, opt->rectype[0], isb_by_prn)) {
                for (f=0;f<nf;f++) meas[f+nf]=meas[f+nf]-isb_by_prn[f]; /* P */
             } else {
                for (f=0;f<nf;f++) meas[f+nf]=meas[f+nf]-isb[f][1]; /* P */
            }
        }
        pt0[sat-1]=nav->ssr[sat-1].t0[5];


#ifdef ENA_PPP_RTK
        /* 1/4 cycle phase shift correction in case of L2C */
        if (opt->phasshft==ISBOPT_TABLE&&isL2C(obs[i].code[1])) {
            for (j=0;j<nav->sfts.n;j++) {
                if (strcmp(nav->sfts.rectyp[j],opt->rectype[0])==0) {
                    ydif=nav->sfts.bias[j];
                }
            }
            meas[1]=meas[1]+ydif*lam[1];
        }
        
        /* measurementd phase and code */
        for (f=0;f<nf;f++) {
            if (lam[f]==0.0||obs[i].L[f]==0.0||obs[i].P[f]==0.0) continue;
            if (obs[i].SNR[f]>0.0&&testsnr(0,f,azel[i*2+1],obs[i].SNR[f]*0.25,&opt->snrmask)) {
                trace(2,"snr error: sat=%2d f=%d el=%.2f SNR=%.2f\n",
                      sat,f,azel[i*2+1]*R2D,obs[i].SNR[f]*0.25);
                continue;
            }
            meas[f   ]+=obs[i].L[f]*lam[f];
            meas[f+nf]+=obs[i].P[f];
        }
        
        qj=selfreqpair(sat,opt,obs+i);
        for (j=0;j<2;j++) { /* for phase and code */
            for (f=0;f<nf;f++) {
                if (meas[nf*j+f]==0.0) continue;
                if (j==0) ssat[sat-1].code[f]=obs[i].code[f];
                if (f>0&&!(f&qj)) continue;
                if (f==1&&(satsys(sat,NULL)==SYS_GAL)) continue;
                if (j==0&&(osr[i].pbias[f]==CSSRINVALID)) {
                    trace(2,"invalid pbias sat=%2d f=%1d\n",sat,f);
                    continue;
                }
                if (j==1&&(osr[i].cbias[f]==CSSRINVALID)) {
                    trace(2,"invalid cbias sat=%2d f=%1d \n",sat,f);
                    continue;
                }
                y[nf*i*2+nf*j+f]=meas[nf*j+f]-modl[nf*j+f];
            }
            nv+=nf;
        }
#endif
      
#ifndef CSSR2OSR_VRS
        if (sys == SYS_GAL) {
            osr[i].antr[1] = 0.0;
            osr[i].wupL[1] = 0.0;
            osr[i].CPC[1]  = 0.0;
            osr[i].PRC[1]  = 0.0;
        }
#else
        obs[k].time=obs[i].time;
        obs[k].sat=sat;
        
        for (j=0;j<nsig;j++) {
            obs[k].code[j] = nav->ssr[sat-1].smode[j];
            obs[k].SNR[j] = 40.0/0.25;
        }
        /* clear observations that are not available */
        for (j=nsig;j<NFREQ+NEXOBS;j++) {
            obs[k].P[j] = 0;
            obs[k].L[j] = 0;
            obs[k].LLI[j] =0;
            obs[k].code[j] = 0;
            obs[k].SNR[j] = 0;
        }
        obs[k].osr_i=i;
        k++;
#endif
        if (osrlog) {
            tow = time2gpst(obs[i].time, NULL);
            getorbitclock(tow, sat, &osr[i].orb, &osr[i].clk);
            trace(3, "OSRRES,%.1f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"
                  "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"
                  "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.9f,%.9f,%.3f\n",
                  tow, sys, prn,
                  osr[i].pbias[0], osr[i].pbias[1], osr[i].pbias[2],
                  osr[i].cbias[0], osr[i].cbias[1], osr[i].cbias[2],
                  osr[i].trop,
                  osr[i].iono,
                  osr[i].antr[0], osr[i].antr[1], osr[i].antr[2],
                  osr[i].relatv,
                  osr[i].wupL[0], osr[i].wupL[1], osr[i].wupL[2],
                  osr[i].compL[0], osr[i].compL[1], osr[i].compL[2],
                  osr[i].sis,
                  osr[i].CPC[0], osr[i].CPC[1], osr[i].CPC[2],
                  osr[i].PRC[0], osr[i].PRC[1], osr[i].PRC[2],
                  osr[i].orb, osr[i].clk,
                  pos[0]*R2D,pos[1]*R2D,pos[2]);
        }
#ifdef ENA_PPP_RTK
        if (osrlog && fp_osr) {
            tow = time2gpst(obs[i].time, NULL);
            getorbitclock(tow, sat, &osr[i].orb, &osr[i].clk);
            fprintf(fp_osr, "OSRRES,%.1f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"
                    "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"
                    "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.9f,%.9f,%.3f\n",
                    tow, sys, prn,
                    osr[i].pbias[0], osr[i].pbias[1], osr[i].pbias[2],
                    osr[i].cbias[0], osr[i].cbias[1], osr[i].cbias[2],
                    osr[i].trop,
                    osr[i].iono,
                    osr[i].antr[0], osr[i].antr[1], osr[i].antr[2],
                    osr[i].relatv,
                    osr[i].wupL[0], osr[i].wupL[1], osr[i].wupL[2],
                    osr[i].compL[0], osr[i].compL[1], osr[i].compL[2],
                    osr[i].sis,
                    osr[i].CPC[0], osr[i].CPC[1], osr[i].CPC[2],
                    osr[i].PRC[0], osr[i].PRC[1], osr[i].PRC[2],
                    osr[i].orb, osr[i].clk,
                    pos[0]*R2D,pos[1]*R2D,pos[2]);
        }
#endif
    }

#ifdef ENA_PPP_RTK
    return nv;
#else
    return k;
#endif
}

/* compensate time variation of signal bias and ionosphere delay  */
extern void compensatedisp(const nav_t *nav,const int *index,
                           const obsd_t *obs, int sat,
                           const double iono, const double *pb,
                           double *compL, int *pbreset, const prcopt_t *opt,
                           const ssat_t ssat)
{
    gtime_t time=obs->time;
    static gtime_t t0[MAXSAT]={0},tm[MAXSAT]={0};
    int i,k,qi,qj,isat=sat-1,oft=isat*NFREQ;
    static double b0[MAXSAT*NFREQ]={0},bm[MAXSAT*NFREQ]={0};
    static double iono0[MAXSAT]={0},ionom[MAXSAT]={0},coef[MAXSAT]={0};
    static int slip[MAXSAT*NFREQ]={0};
    const double *lam=nav->lam[obs->sat-1];
    double disp0,dispm,dt,dgf,fi;
    int nf=opt->nf,flag;
    int nsig=nav->ssr[sat-1].nsig;
    
#ifdef CSSR2OSR_VRS
    nf=nsig;
#endif

    /* long interval */
    for (k=flag=0;k<nav->stec[index[0]].n;k++) {
        if (sat==nav->stec[index[0]].data[k].sat) {
            flag=1;
            break;
        }
    }
    if (flag==1&&timediff(nav->stec[index[0]].data[k].time,t0[isat])>0.0) {
        if (opt->posopt[5]==1) {
            tm[isat]=t0[isat];
            t0[isat]=nav->stec[index[0]].data[k].time;
            dt=timediff(t0[isat],tm[isat]);
            if (dt<=0.0) return;

            for (i=0;i<nf;i++) {
                if (b0[oft+i]==0.0||iono0[isat]==0.0) {
                    b0[oft+i]=pb[i];iono0[isat]=iono;
                    return;
                }
            }

            /* signal bias */
            for (i=0;i<nf;i++) {
                bm[oft+i]=b0[oft+i];
                b0[oft+i]=pb[i];
            }
            /* ionosphere */
            ionom[isat]=iono0[isat];
            iono0[isat]=iono;

            for (i=1;i<nf;i++) {
                qi=0;qj=i;fi=lam[qj]/lam[qi];
                if (pb[qi]==CSSRINVALID||pb[qj]==CSSRINVALID||iono==0.0) continue;
                dispm=-FREQ2/FREQ1*(1.0-fi*fi)*ionom[isat]+bm[oft+qi]-bm[oft+qj];
                disp0=-FREQ2/FREQ1*(1.0-fi*fi)*iono0[isat]+b0[oft+qi]-b0[oft+qj];
                coef[isat+(i-1)*MAXSAT]=(disp0-dispm)/dt;
            }
        } else {
            for (i=0;i<nf;i++) {
                b0[oft+i]=obs->L[i]*lam[i];
                slip[oft+i]=0;
            }
            tm[isat]=t0[isat];
            t0[isat]=nav->stec[index[0]].data[k].time;
        }
    }

    dt=timediff(time,t0[isat]);
    if (opt->posopt[5]==1) {
        for  (i=1;i<nf;i++) {
            qi=0;qj=i;fi=lam[qj]/lam[qi];
            if (fabs(coef[isat+(i-1)*MAXSAT]/(FREQ2/FREQ1*(1.0-fi*fi)))>0.008) return;
            if (pbreset[qi]||pbreset[qj]) {
                coef[isat]=0.0;
                return;
            }
            compL[qi]=compL[qi]==0.0?(1.0/(1.0-SQR(fi))*coef[isat+(i-1)*MAXSAT]*dt):compL[qi];
            compL[qj]=SQR(fi)/(1.0-SQR(fi))*coef[isat+(i-1)*MAXSAT]*dt;
        }
    } else {
        for (i=0;i<nf;i++) {
            if (ssat.slip[i]>0) slip[oft+i]=1;
        }
        for  (i=1;i<nf;i++) {
            qi=0;qj=i;fi=lam[qj]/lam[qi];
            if (slip[oft+qi]||slip[oft+qj]||
                    obs->L[qi]*lam[qi]==0.0||obs->L[qj]*lam[qj]==0.0||pbreset[qi]||pbreset[qj]) continue;
            dgf=obs->L[qi]*lam[qi]-obs->L[qj]*lam[qj]-(b0[oft+qi]-b0[oft+qj]);
            compL[qi]=compL[qi]==0.0?(1.0/(1.0-SQR(fi))*dgf):compL[qi];
            compL[qj]=SQR(fi)/(1.0-SQR(fi))*dgf;
        }
    }
}

/* ionosphere and antenna corrected measurements -----------------------------*/
extern int corrmeas(const obsd_t *obs, nav_t *nav, const double *pos,
                    const double *azel, const prcopt_t *opt,
                    const int *index, const int n, const double *weight,
                    const double *Gmat, const double *Emat, const ssat_t ssat,
                    int *brk, osrd_t *osr, int *pbreset, int sisadjust)
{
    const double *lam=nav->lam[obs->sat-1];
    double vari,dant[NFREQ+NEXOBS]={0},compL[NFREQ+NEXOBS]={0};
    double stec=0.0,quality=0.0, rate, t5, t6;
    double pbias[NFREQ+NEXOBS]={0},cbias[NFREQ+NEXOBS]={0};
    int i,j,sat,smode,nsig,nf=opt->nf;
    int flag;
    double dt;
    static gtime_t currtime = {0};

    trace(3,"corrmeas: time=%s, sat=%2d\n",time_str(obs->time,0),obs->sat);

    sat=obs->sat;
    for (i=0;i<nf;i++) pbias[i]=cbias[i]=CSSRINVALID;
    
    /* decode phase and code bias */
    t5=timediff(obs->time,nav->ssr[sat-1].t0[4]); /* cbais */
    t6=timediff(obs->time,nav->ssr[sat-1].t0[5]); /* pbais */

    if (fabs(t5)>MAXAGESSR_BIAS||fabs(t6)>MAXAGESSR_BIAS) {
        trace(2,"age of ssr bias error: time=%s sat=%2d tc(cbias),tp(pbias)=%.0f,%.0f\n",
            time_str(obs->time,0),sat,t5,t6);
        memset(&nav->ssr[sat-1].t0[4], 0x00, sizeof(nav->ssr[sat-1].t0[4]));
        memset(&nav->ssr[sat-1].t0[5], 0x00, sizeof(nav->ssr[sat-1].t0[5]));
        return 0;
    }

    /* decode phase and code bias */
#ifndef CSSR2OSR_VRS
    switch (nav->rtcmmode) {
        case RTCMMODE_CSSR:
            nsig=nav->ssr[sat-1].nsig;
            for (i=0;i<nf;i++) {
                for (j=0;j<nsig;j++) {
                    smode=nav->ssr[sat-1].smode[j];
                    if (obs->code[i]==smode) {
                        pbias[i]=nav->ssr[sat-1].pbias[smode-1];
                        cbias[i]=nav->ssr[sat-1].cbias[smode-1];
                        break;
                    }
                }
            }
            break;
        default:
            trace(2,"not supported ssr format\n");
            return 0;
    }
#else
    nsig=nav->ssr[sat-1].nsig;   
    for (i=0;i<nsig;i++) {
        pbias[i]=cbias[i]=CSSRINVALID;
        smode=nav->ssr[sat-1].smode[i];
        pbias[i]=nav->ssr[sat-1].pbias[smode-1];
        cbias[i]=nav->ssr[sat-1].cbias[smode-1];
    }
#endif

    /* ionosphere correction */
    if (!stec_grid_data(nav,index,obs->time,sat,n,weight,Gmat,Emat,&stec,&rate,&vari,&quality,brk)) {
        trace(2,"corrmea: iono correction error: time=%s sat=%2d ionoopt=%d\n",
              time_str(obs->time,2),obs->sat,opt->ionoopt);
        return 0;
    }
    
    /* inconsistency check */
    for (i=flag=0;i<nav->stec[index[0]].n;i++) {
        if (sat==nav->stec[index[0]].data[i].sat) {
            flag=1;
            break;
        }
    }
    if (flag==1) {
        dt=timediff(nav->stec[index[0]].data[i].time,nav->ssr[sat-1].t0[5]);
        nav->ssr[sat-1].t0[8]=nav->stec[index[0]].data[i].time;
        if (dt<0.0||dt>=30.0) {
            trace(2,"inconsist ssr correction (iono-pbias):time=%s sat=%2d dt=%.2f\n",
                  time_str(obs->time,0),sat,dt);
            return 0;
        }
        dt=timediff(nav->stec[index[0]].data[i].time,nav->ssr[sat-1].t0[0]);
        if (dt<-30.0||dt>30.0) {
            trace(2,"inconsist ssr correction (iono-orbit):time=%s sat=%2d dt=%.2f\n",
                  time_str(obs->time,0),sat,dt);
            return 0;
        }
    }
    dt=timediff(nav->ssr[sat-1].t0[5],nav->ssr[sat-1].t0[4]);
    trace(5, "corrmeas(): pbias_tow=%.1f, cbias_tow=%.1f, dt=%.1f\n", time2gpst(nav->ssr[sat-1].t0[5], NULL), time2gpst(nav->ssr[sat-1].t0[4], NULL), dt);
    if (opt->posopt[11] != 0) {
        if ((!opt->posopt[9]&&fabs(dt)>0.0)||(opt->posopt[9]&&(dt>0.0||dt<-30.0))) {
            trace(2,"inconsist correction (pbias-cbias):time=%s sat=%2d dt=%.2f\n",
                  time_str(obs->time,0),sat,dt);
            return 0;
        }
    } else {
        if ((sisadjust&&fabs(dt)>0.0)||(!sisadjust&&(dt>0.0||dt<-30.0))) {
            trace(2,"inconsist correction (pbias-cbias):time=%s sat=%2d dt=%.2f\n",
                  time_str(obs->time,0),sat,dt);
            return 0;
        }
    }
    
    if (timediff(obs->time, currtime) != 0.0) {
        trace(2, "data time, obs, %.1f, orbit, %.1f, clock, %.1f, pbias, %.1f, cbias, %.1f, iono, %.1f\n", time2gpst(obs->time, NULL), time2gpst(nav->ssr[sat-1].t0[0], NULL),
              time2gpst(nav->ssr[sat-1].t0[1], NULL), time2gpst(nav->ssr[sat-1].t0[5], NULL), time2gpst(nav->ssr[sat-1].t0[4], NULL),
              time2gpst(nav->stec[index[0]].data[i].time, NULL));
        currtime = obs->time;
    }

    /* compensate time variation of signal bias and ionosphere delay */
    if (opt->posopt[5]>0) {
        compensatedisp(nav,index,obs,sat,stec,pbias,compL,pbreset,opt,ssat);
    }

    /* receiver antenna phase center correction */
    antmodel(opt->pcvr,opt->antdel[0],azel,opt->posopt[1],dant);

    /* ionosphere and windup corrected phase and code */

     /* for (i=0;i<nf;i++) { ssr2osr用コード */
#ifdef CSSR2OSR_VRS
    for (i=0;i<nsig;i++){
#else
    for (i=0;i<nf;i++) {
#endif
        osr->cbias[i]=cbias[i];
        osr->pbias[i]=pbias[i];
        osr->compL[i]=compL[i];
    }

    for (i=0;i<nf;i++) {
       osr->wupL[i]=lam[i]*ssat.phw;
       osr->antr[i]=dant[i];
    }
    osr->iono=stec;
    
    return 1;
}

/* initialize PPP-RTK control -------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/

extern void rtkinitppprtk(rtk_t *rtk, const prcopt_t *opt)
{
    sol_t sol0={{0}};
    ambc_t ambc0={{{0}}};
    ssat_t ssat0={0};
    int i;

    trace(3,"rtkinitppprtk :\n");

    rtk->sol=sol0;
    rtk->nx=NX(opt);
    rtk->x=zeros(rtk->nx,1);
    rtk->nfix=rtk->neb=0;
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i]=ssat0;
    }
    for (i=0;i<MAXERRMSG;i++) rtk->errbuf[i]=0;
    rtk->opt=*opt;
    rtk->sisadjust = 1;

#ifndef CSSR2OSR_VRS
    for (i=0;i<6;i++) rtk->rb[i]=0.0;
    rtk->na=NR(opt);
    rtk->tt=0.0;
    rtk->P=zeros(rtk->nx,rtk->nx);
    rtk->Q=zeros(rtk->nx,rtk->nx);
    rtk->xa=zeros(rtk->na,1);
    rtk->Pa=zeros(rtk->na,rtk->na);
    for (i=0;i<MAXSAT;i++) {
        rtk->ambc[i]=ambc0;
    }
#endif
}

extern void rtkfreeppprtk(rtk_t *rtk)
{
    trace(3,"rtkfreeppprtk :\n");

    rtk->nx = 0;
    free(rtk->x); rtk->x = NULL;
#ifndef CSSR2OSR_VRS
    rtk->na = 0;
    free(rtk->P); rtk->P = NULL;
    free(rtk->Q); rtk->Q = NULL;
    free(rtk->xa); rtk->xa = NULL;
    free(rtk->Pa); rtk->Pa = NULL;
#endif
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
extern int ssr2osr(rtk_t *rtk, const obsd_t *obs, const int n, nav_t *nav,
                 osrd_t *osr, const int mode)

    {
    double *rs,*dts,*var,*e,*azel, *cpc_=NULL;
    gtime_t *pt0_=NULL;
    double dt,dt_p;
    int i,f,k=0,sat,j,nf=rtk->opt.nf,sati,ephopt;
    int svh[MAXOBS];
    float age;
    static grid_t grid;
    static int backup = FALSE;
    double pos[3];
    gtime_t temp;
    int nn;

    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);

    e = mat(n,3);

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys = satsys(i+1,NULL);
        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat[j] = rtk->ssat[i].snr[j] = 0;
    }

    for (i=0;i<MAXSAT;i++) rtk->ssat[i].fix[0]=0;

#ifdef CSSR2OSR_VRS
    for (i=0;i<3;i++) rtk->x[i]=rtk->opt.ru[i];
#else
    /* temporal update of states */
    if (rtk->opt.mode==PMODE_SSR2OSR_FIXED) {
        for (i=0;i<3;i++) rtk->x[i]=rtk->opt.ru[i];
    } else {
        if (norm(rtk->sol.rr,3)>0.0) {
            for (i=0;i<3;i++) rtk->x[i]=rtk->sol.rr[i];
        }
    }
#endif

    grid.network = rtk->opt.netnum;
    ecef2pos(rtk->x, pos);        
    if ((nn = get_grid_index(nav, pos, &grid, &(rtk->opt), obs[0].time, FALSE)) == 0) {
        if (backup == FALSE || timediff(obs[0].time, (temp = GetBackupCSSRTime())) > 180.0) {
            free(rs); free(dts); free(var); free(e); free(azel);
            rtk->sol.stat = SOLQ_SINGLE;
            return 0;
        } else {
            trace(3, "use backup CSSR data: obs=%.1f, CSSR=%.1f, age=%.1f\n",
                time2gpst(obs[0].time, NULL), time2gpst(temp, NULL),
                timediff(obs[0].time, temp));
            RestoreCurrentCSSR(&grid);
        }
    } else {
        if (GetCloseCSSR(obs[0].time, grid.network) == FALSE) {
            free(rs); free(dts); free(var); free(e); free(azel);
            rtk->sol.stat = SOLQ_SINGLE;
            return 0;
        }
        BackupCurrentCSSR(&grid);
        backup = TRUE;
    }
    
    nav->facility = GetCurrentCSSRFacility();
    nav->ssrtime = GetCurrentCSSRTime();
    CheckCSSRFacility(nav, grid.network);
    
    if (rtk->opt.posopt[11] != 0) {
        if (!rtk->opt.posopt[9] && IsSISAdjust()) {
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

    /* ssr age */
    rtk->sol.age=1e4;
    for (i=0;i<n&&i<MAXOBS;i++) {
        sati = obs[i].sat;
        age=timediff(obs[i].time,nav->ssr[sati-1].t0[1]);
        if (rtk->sol.age>age) rtk->sol.age=age;
    }

    /* satellite positions and clocks */
    saveposition(rtk->x);
    satposs(obs[0].time,obs,n,nav,EPHOPT_SSRAPC,rs,dts,var,svh);
    for (i=0;i<MAXSAT;i++) {
        satpos_ssr_sis(obs[0].time,obs[0].time,i+1,rtk,nav);
    }

    /* undifferential residuals */
    if (!(k=zdres(obs,n,rs,dts,var,svh,nav,rtk->x,NULL,e,azel,rtk,TRUE,cpc_,pt0_,&grid,rtk->ssat,&rtk->opt,&rtk->sol,osr))) {
        trace(2,"rover initial position error\n");
    }

    if ((nav->filreset == TRUE) && (rtk->opt.mode==PMODE_SSR2OSR))  nav->filreset = FALSE;

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
            if (rtk->ssat[i].fix[j]==2) rtk->ssat[i].fix[j]=1;
            if (rtk->ssat[i].slip[j]&1) rtk->ssat[i].slipc[j]=1;
        }
    }

    free(rs); free(dts); free(var); free(e); free(azel);
    return k;

}
