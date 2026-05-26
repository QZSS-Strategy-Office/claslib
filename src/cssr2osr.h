/*------------------------------------------------------------------------------
* cssr2osr.h : Compact SSR constants, types and function prototypes
*
*          Copyright (C) 2015- by Mitsubishi Electric Corporation, All rights reserved.
*-----------------------------------------------------------------------------*/
#ifndef CSSR2OSR_H
#define CSSR2OSR_H

#define SQR(x)      ((x)*(x))
#define MAX_NGRID 4                 /* number of gridded points for interpolation */
#define MAXPBCORSSR 20.0            /* max phase bias correction of ssr (m) */
#define CSSRINVALID -10000          /* invalid value */

#define NP(opt)     ((opt)->dynamics==0?3:9) /* number of pos solution */
#define NF(opt)     ((opt)->nf)
#define IC(s,opt)   (NP(opt)+(s))      /* state index of clocks (s=0:gps,1:glo) */
#define IT(opt)     (IC(0,opt)+NSYS)   /* state index of tropos */
#define NI(opt)     ((((opt)->ionoopt!=IONOOPT_EST)&&((opt)->ionoopt!=IONOOPT_EST_ADPT))?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?1:3))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))  /* number of estimated states */
#define II(s,opt)   (NP(opt)+(s)-1)    /* ionos (s:satellite no) */
#define ITT(opt)    (NP(opt)+NI(opt))  /* tropos (r:0=rov,1:ref) */
#define IB(s,f,opt)   (NR(opt)+MAXSAT*(f)+(s)-1)    /* state index of phase bias */

extern double prectrop(gtime_t time, const double *pos, const double *azel,
                       const prcopt_t *opt, const double zwd, double ztd);
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
                 ssat_t *ssat, prcopt_t *opt, sol_t *sol, osrd_t *osr);
extern int CheckSelectGridData(nav_t *nav, gtime_t obstime, int network, int num, int *index);
extern void compensatedisp(const nav_t *nav,const int *index,
                           const obsd_t *obs, int sat,
                           const double iono, const double *pb,
                           double *compL, int *pbreset, const prcopt_t *opt,
                           const ssat_t ssat);
extern int corrmeas(const obsd_t *obs, nav_t *nav, const double *pos,
                    const double *azel, const prcopt_t *opt,
                    const int *index, const int n, const double *weight,
                    const double *Gmat, const double *Emat, const ssat_t ssat,
                    int *brk, osrd_t *osr, int *pbreset, int sisadjust);
extern void rtkinitppprtk(rtk_t *rtk, const prcopt_t *opt);
extern void rtkfreeppprtk(rtk_t *rtk);
extern int ddres(rtk_t *rtk, const nav_t *nav, const double *x,
                 const double *P, const obsd_t *obs, double *y, double *e,
                 double *azel, int n, double *v, double *H, double *R,
                 int *vflg, int niter);
extern void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);
extern int selfreqpair(const int sat, const prcopt_t *opt,const obsd_t *obs);
#endif /* CSSR2OSR_H */
