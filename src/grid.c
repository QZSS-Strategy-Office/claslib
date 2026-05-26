/*
 * grid.c: grid correction
 *
 *          Copyright (C) 2015- by Mitsubishi Electric Corporation, All rights reserved.
 */

#include "rtklib.h"
#include "clasgrid.h"


#define MAX_NGRID   4           /* number of grids for interpolation */
#define MAX_GRID_CASHE   20      /* grids cashe */
#define MAX_DIST    120.0       /* max distance to grid (km) */
#define MAX_AGE     60.0       /* max age of difference (s) */
#define VAR_NOTEC   SQR(30.0)   /* variance of no tec */
#define MIN_EL      0.0         /* min elevation angle (rad) */
#define MIN_HGT     -1000.0     /* min user height (m) */
#define NUM_GRID    21          /* number of gridded points for each area 32*/

#define SQR(x)      ((x)*(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))

#define CSSRINVALID -10000      /* invalid value*/


/* add stec data -------------------------------------------------------------*/
extern int add_data_stec(stec_t *stec, gtime_t time, int sat, int slip,
                         double iono, double rate, double rms, double quality)
{
    stec->data[stec->n].flag=1;
    stec->data[stec->n].time=time;
    stec->data[stec->n].sat=(unsigned char)sat;
    stec->data[stec->n].slip=(unsigned char)slip;
    stec->data[stec->n].iono=(float)iono;
    stec->data[stec->n].rate=(float)rate;
    stec->data[stec->n].quality=(float)quality;
    stec->data[stec->n++].rms=(float)rms;
    return 1;
}
/* add stec data -------------------------------------------------------------*/
extern int add_data_trop(zwd_t *z, gtime_t time, double zwd, double ztd,
                         double quality, double rms, int valid)
{
    z->data[z->n].time=time;
    z->data[z->n].valid=(unsigned char)valid;
    z->data[z->n].zwd=(float)zwd;
    z->data[z->n].ztd=(float)ztd;
    z->data[z->n].quality=(float)quality;
    z->data[z->n++].rms=(float)rms;
    
    return 1;
}
typedef struct _GridInfo {
    double pos[3];
    double weight;
    int network;
    int index;
} GridInfo;
int gridsel = 2;


static int PickupCandidatesOfCSSRGrid(int network, const double *pos, GridInfo *gridinfo)
{
    int start = (network > 0 ? network: 1), end = (network > 0 ? network + 1: CSSR_MAX_NETWORK);
    int inet, i, j, k, n = 0;
    double d, dd[2];

    for (inet = start; inet < end; ++inet) {
        for (i = 0; i < RTCM_SSR_MAX_GP; ++i) {
            if (fabs(-1.0 - clas_grid[inet][i][0]) < 0.0001 && fabs(-1.0 - clas_grid[inet][i][1]) < 0.0001 && fabs(-1.0 - clas_grid[inet][i][2]) < 0.0001) {
                break;
            }
            if (fabs(clas_grid[inet][i][2]) > 0.0001) {
                continue;
            }
            
            /* distance to grid (m) */
            dd[0] = RE_WGS84 * (pos[0] - clas_grid[inet][i][0] * D2R);
            dd[1] = RE_WGS84 * (pos[1] - clas_grid[inet][i][1] * D2R) * cos(pos[0]);
            if ((d = MAX(norm(dd, 2), 1.0)) > MAX_DIST * 1000.0) {
                continue;
            }
            
            trace(5, "get_grid_index: inet=%2d, pos=%3.3f %3.3f, dist=%6.3f[km]\n", inet,
                clas_grid[inet][i][0], clas_grid[inet][i][1], d / 1000.0);
            
            if (n <= 0) {
                gridinfo[n].pos[0] = clas_grid[inet][i][0];
                gridinfo[n].pos[1] = clas_grid[inet][i][1];
                gridinfo[n].pos[2] = clas_grid[inet][i][2];
                gridinfo[n].network = inet;
                gridinfo[n].weight = d;
                gridinfo[n].index = i;
                ++n;
            } else {
                for (j = 0; j < n; j++) if (d < gridinfo[j].weight) break;
                if (j >= MAX_GRID_CASHE) continue;
                for (k = MIN(n, MAX_GRID_CASHE - 1); k > j; k--) {
                    gridinfo[k].pos[0] = gridinfo[k-1].pos[0];
                    gridinfo[k].pos[1] = gridinfo[k-1].pos[1];
                    gridinfo[k].pos[2] = gridinfo[k-1].pos[2];
                    gridinfo[k].network = gridinfo[k-1].network;
                    gridinfo[k].weight = gridinfo[k-1].weight;
                    gridinfo[k].index = gridinfo[k-1].index;
                }
                gridinfo[j].pos[0] = clas_grid[inet][i][0];
                gridinfo[j].pos[1] = clas_grid[inet][i][1];
                gridinfo[j].pos[2] = clas_grid[inet][i][2];
                gridinfo[j].network = inet;
                gridinfo[j].weight = d;
                gridinfo[j].index = i;
                if (n < MAX_GRID_CASHE) n++;
            }
        }
    }
    return n;
}

static int FindGridSurroundPos(GridInfo **gridindex, int n, GridInfo *gridinfo)
{
    int first_index = (gridinfo[1].network == gridinfo[2].network && gridinfo[0].network != gridinfo[1].network ? 1 : 0);
    double dd12[2], dd13[2], dd14[2];
    GridInfo *grid2, *grid3;
    int flag, i;
    
    for (i = 1, gridindex[0] = &gridinfo[first_index], flag = FALSE; i < n; ++i) {
        if (gridinfo[first_index].network == gridinfo[i].network && first_index != i) {
            dd12[0] = gridinfo[i].pos[0]*D2R - gridinfo[first_index].pos[0]*D2R;
            dd12[1] = gridinfo[i].pos[1]*D2R - gridinfo[first_index].pos[1]*D2R;
            gridindex[1] = &gridinfo[i];
            grid2 = &gridinfo[i];
            flag = TRUE;
            break;
        }
    }
    
    if (flag == FALSE) {
        trace(2, "can't find second grid: network=%d\n", gridinfo[first_index].network);
        return 1;
    }
    
    for (i = 1, flag = FALSE; i < n; i++) {
        if (gridindex[1]->index == gridinfo[i].index || gridinfo[first_index].network != gridinfo[i].network) {
            continue;
        }
        dd13[0] = gridinfo[i].pos[0]*D2R - gridinfo[first_index].pos[0]*D2R;
        dd13[1] = gridinfo[i].pos[1]*D2R - gridinfo[first_index].pos[1]*D2R;
        if (fabs(dot(dd12, dd13, 2) / (norm(dd12, 2) * norm(dd13, 2))) < 1.0) {
            gridindex[2] = &gridinfo[i];
            grid3 = &gridinfo[i];
            flag = TRUE;
            break;
        }
    }
    
    if (flag == FALSE) {
        trace(2, "can't find third grid: network=%d\n", gridinfo[first_index].network);
        return 1;
    }
    
    for (i = 1, flag = FALSE; i < n; i++) {
        if (gridindex[1]->index == gridinfo[i].index || gridindex[2]->index == gridinfo[i].index || gridinfo[first_index].network != gridinfo[i].network) {
            continue;
        }
        dd14[0] = gridinfo[i].pos[0]*D2R - gridinfo[first_index].pos[0]*D2R;
        dd14[1] = gridinfo[i].pos[1]*D2R - gridinfo[first_index].pos[1]*D2R;
        if (fabs(dot(dd12, dd14, 2) / (norm(dd12, 2) * norm(dd14, 2))) < 1.0 &&
            fabs(dot(dd13, dd14, 2) / (norm(dd13, 2) * norm(dd14, 2))) < 1.0) {
            if ((fabs(grid2->pos[1] - gridinfo[i].pos[1]) < 0.04 && fabs(grid3->pos[0] - gridinfo[i].pos[0]) < 0.04) ||
                (fabs(grid2->pos[0] - gridinfo[i].pos[0]) < 0.04 && fabs(grid3->pos[1] - gridinfo[i].pos[1]) < 0.04)) {
                gridindex[3] = &gridinfo[i];
                flag = TRUE;
                break;
            }
        }
    }
    
    if (flag == FALSE) {
        trace(2, "can't find fourth grid: network=%d\n", gridinfo[first_index].network);
        n = 3;
    } else {
        n = 4;
    }
    return n;
}

static void Calculation3PointsWeight(double *weight, const double *pos, GridInfo **gridindex)
{
    double *U, *E, dd[2];
    int i;
    
    U=mat(2,2); E=mat(2,1);
    for (i=1;i<3;i++) {
        U[(i-1)*2+0] = gridindex[i]->pos[0]-gridindex[0]->pos[0];
        U[(i-1)*2+1] = gridindex[i]->pos[1]-gridindex[0]->pos[1];
    }
    trace(5,"Umat=\n"); tracemat(5,U,2,2,10,5);
    
    /* UD decomposition */
    if (matinv(U,2)) trace(2,"calculation error of G inverse\n");
    trace(5,"Umat=\n"); tracemat(5,U,2,2,10,5);
    
    E[0] = pos[0]*R2D-gridindex[0]->pos[0];
    E[1] = pos[1]*R2D-gridindex[0]->pos[1];
    matmul("NN",2,1,2,1.0,U,E,0.0,dd);
    free(U); free(E);
    
    weight[0] = 1.0 - dd[0] - dd[1];
    weight[1] = dd[0];
    weight[2] = dd[1];
}

static int IsExistsInside(const double *pos, int n, GridInfo **gridindex, int rtcmmode, GridInfo *nearestgrid)
{
    int flag = 0xff;
    
    if (flag == 0x00) {
        if (rtcmmode == RTCMMODE_CSSR) {
            if (nearestgrid->weight < gridindex[0]->weight) {
                trace(2, "change to the nearest grid: inet=%d, index=%d, pos=%6.3f %6.3f %6.3f, distance=%.2fKm\n",
                    nearestgrid->network, nearestgrid->index, nearestgrid->pos[0], nearestgrid->pos[1],
                    nearestgrid->pos[2], nearestgrid->weight / 1000.0);
                gridindex[0] = nearestgrid;
                return 1;
            }
        }
        trace(2, "change to the nearest grid: inet=%d, index=%d, pos=%6.3f %6.3f %6.3f, distance=%.2fKm\n",
            gridindex[0]->network, gridindex[0]->index, gridindex[0]->pos[0], gridindex[0]->pos[1],
            gridindex[0]->pos[2], gridindex[0]->weight / 1000.0);
        n = 1;
    }
    return n;
}

static int CheckGridStatus(gtime_t obstime, int n, GridInfo **gridindex, int rtcmmode)
{
    int valid = (n == 4 ? 0x0f: (n == 3 ? 0x07: 0x01));
    int flag, change=0, i, j;
    GridInfo *temp[4];
    
    for (i = j = 0, flag = 0x00; i < n; ++i) {
        if (CheckGridData(obstime, gridindex[i]->network, gridindex[i]->index) == TRUE) {
            temp[j++] = gridindex[i];
            flag |= (1 << i);
        }
    }
    
    if (flag != valid && flag > 0x00) {
        switch (rtcmmode) {
        case RTCMMODE_CSSR:
            change = (timediff(obstime, GetBackupCSSRTime()) < SSRVALIDAGE ? FALSE: TRUE);
            break;
        default:
            break;
        }
        if (change && temp[0]->weight < 60000.0) {    /* need option ? */
            trace(2, "change to the nearest grid: inet=%d, index=%d, pos=%6.3f %6.3f %6.3f, distance=%.2fKm\n",
                temp[0]->network, temp[0]->index, temp[0]->pos[0], temp[0]->pos[1], temp[0]->pos[2],
                temp[0]->weight / 1000.0);
                gridindex[0] = temp[0];
                return 1;
        } else {
            return 0;
        }
    } else if (flag == 0x00) {
        return 0;
    }
    return n;
}

static void CalculationGridWeight(double *weight, const double *pos, int n, GridInfo **gridindex)
{
    double sum;
    int i;
    
    if (n != 3) {
        sum=0.0;
        for (i=0;i<n;i++) sum += 1.0/gridindex[i]->weight;
        for (i=0;i<n;i++) weight[i] = 1.0/(gridindex[i]->weight*sum);
    } else {
        Calculation3PointsWeight(weight, pos, gridindex);
    }
}

static void OutputSelectedGrid(const double *pos, int n, GridInfo **gridindex, double *weight, const nav_t *nav)
{
    static int savenetwork[4];
    static int saveindex[4];
    static int savenum = -1;
    int flag = FALSE;
    int i;
    
    for (i = 0; i < savenum; ++i) {
        if (savenetwork[i] != gridindex[i]->network || saveindex[i] != gridindex[i]->index) {
            flag = TRUE;
            break;
        }
    }
    if (savenum == n && flag == FALSE) {
        return;
    }
    
    trace(1, "pos est: pos=%6.3f %6.3f               \n", pos[0]/D2R, pos[1]/D2R);
    for (i = 0; i < n; i++) {
        trace(1, "selected_grid%d: inet=%d, weight=%.2f, pos=%6.3f %6.3f, dist=%.2fkm, index=%d\n",
            i + 1, gridindex[i]->network, weight[i], gridindex[i]->pos[0], gridindex[i]->pos[1],
            gridindex[i]->weight / 1000.0, gridindex[i]->index);
    }
    
    for (i = 0; i < n; ++i) {
        savenetwork[i] = gridindex[i]->network;
        saveindex[i] = gridindex[i]->index;
    }
    savenum = n;
}

/* search stec grid -----------------------------------------------------------
* 
* args   : nav_t *nav
*          double *pos
*          int nmax
*          int ofst
*          int num_grid
*          int *index
*          double *weight
*          double *Gmat
*          double *Emat
*          int *netnum
*          prcopt_t *opt
* return : 
* notes  :
*-----------------------------------------------------------------------------*/
extern int get_grid_index(const nav_t *nav, const double *pos, grid_t *grid, prcopt_t *opt, gtime_t obstime, int selectflag)
{
    static GridInfo gridinfo[MAX_GRID_CASHE];
    static GridInfo *gridindex[4];
    double dlat, dlon;
    int i, n = 0;

    trace(3, "get_grid: pos=%.3f %.3f n=%d\n", pos[0]*R2D, pos[1]*R2D, nav->nn);

    if (gridsel >= 2) {
        switch (nav->rtcmmode) {
        case RTCMMODE_CSSR:
            n = PickupCandidatesOfCSSRGrid(grid->network, pos, gridinfo);
            break;
        default:
            break;
        }
    }
    
    if (n >= 3 && gridinfo[0].network <= 12) {
        if (opt->gridsel == 0 || gridinfo[0].weight > (double)opt->gridsel) {
            n = FindGridSurroundPos(&gridindex[0], n, gridinfo);
        } else {
            gridindex[0] = &gridinfo[0];
            n = 1;
        }
    } else if (n == 0) {
        trace(2, "can't find nearby grid: rover pos=%3.3f %3.3f\n",
            pos[0]*R2D, pos[1]*R2D);
        return 0;
    } else {
        gridindex[0] = &gridinfo[0];
        n = 1;
    }

    n = IsExistsInside(pos, n, &gridindex[0], nav->rtcmmode, &gridinfo[0]);
    if (!(n = CheckGridStatus(obstime, n, &gridindex[0], nav->rtcmmode))) {
        trace(2, "can't select valid grid: tow=%.1f\n", time2gpst(obstime, NULL));
        return 0;
    }

    CalculationGridWeight(grid->weight, pos, n, &gridindex[0]);
    OutputSelectedGrid(pos, n, gridindex, grid->weight, nav);

    if (n == 4) {
        double *U=mat(n,n);
        for (i=0;i<n;i++) {
            dlat=gridindex[i]->pos[0]-gridindex[0]->pos[0];
            dlon=gridindex[i]->pos[1]-gridindex[0]->pos[1];
            U[0*n+i]=dlat;U[1*n+i]=dlon;
            U[2*n+i]=dlat*dlon;U[3*n+i]=1.0;
        }
        trace(5,"Umat=\n"); tracemat(5,U,4,n,10,5);

        /* UD decomposition */
        if (matinv(U,n)) trace(2,"calculation error of G inverse\n");
        matcpy(grid->Gmat,U,n,n);
        trace(5,"Gmat=\n"); tracemat(5,grid->Gmat,4,n,10,5);
        free(U);
    }
    dlat=pos[0]*R2D-gridindex[0]->pos[0];
    dlon=pos[1]*R2D-gridindex[0]->pos[1];
    grid->Emat[0]=dlat;grid->Emat[1]=dlon;grid->Emat[2]=dlat*dlon;grid->Emat[3]=1.0;

    grid->network = (grid->network == 0 ? gridindex[0]->network: grid->network);
    for (i = 0; i < n; i++) {
        grid->index[i] = gridindex[i]->index;
    }
    grid->num = n;
    return n;
}

/* search trop data ----------------------------------------------------------*/
extern int trop_data(zwd_t *z, gtime_t time, double *ztd, double *zwd,
                     double *quality, const int idx, int *valid)
{
    double tt;
    int k=0;
    double tdvd,twvd,gh,pos[2];

    trace(4,"trop_data: %s\n",time_str(time,0));

    if (z->n<=0) return 0;

    tt=timediff(time,z->data[k].time);
    if (fabs(tt)>MAXAGESSR_TROP) {
        trace(2,"age of ssr trop error %s tt=%.0f\n",time_str(time,0),tt);
        return 0;
    }
    
    pos[0]=z->pos[0]*D2R;
    pos[1]=z->pos[1]*D2R;
    gh=geoidh(pos);
    if (get_stTv(time,pos[0],0.0,gh,&tdvd,&twvd)) return 0;
    if(z->data[k].ztd==CSSRINVALID||z->data[k].zwd==CSSRINVALID) {
        trace(2,"invalid trop correction %s ztd=%.2f zwd=%.2f\n",
              time_str(time,0),z->data[k].ztd,z->data[k].zwd);
        return 0;
    }
    *ztd=(z->data[k].ztd-z->data[k].zwd)/tdvd;
    *zwd=z->data[k].zwd/twvd;
    *quality=z->data[k].quality;
    *valid=z->data[k].valid;

    return 1;
}
/* ionosphere model by stec grid data ------------------------------------------
* compute ionospheric delay by stec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *iono     O   ionospheric delay (L1) (m)
*          double *rate     O   ionospheric rate (L1) (m/s)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
*          int    brk       O   break flag
* return : status (1:ok,0:error)
* notes  : non-thread-safe
*-----------------------------------------------------------------------------*/
extern int stec_grid_data(const nav_t *nav, const int *index, gtime_t time,
                          int sat, int n, const double *weight,
                          const double *Gmat, const double *Emat, double *iono,
                          double *rate, double *var, double *quality, int *brk)
{
    int i,slip;
    double *ionos,*rates,*rms,*quals;
    double *ionos_,*rates_,*sqrms,*sqrms_,*quals_;

    if (n<=0) return 0;

    ionos=mat(n,1);rates=mat(n,1);rms=mat(n,1);quals=mat(n,1);
    ionos_=mat(n,1);rates_=mat(n,1);sqrms=mat(n,1);sqrms_=mat(n,1);quals_=mat(n,1);

    if (n == 1) {
        if (!stec_data(nav->stec+index[0],time,sat,iono,rate,rms,quals,&slip)) {
            free(ionos);free(rates);free(rms);free(quals);
            free(ionos_);free(rates_);free(sqrms);free(sqrms_);free(quals_);
            return 0;
        }
        if (slip) *brk=1;
        *var=SQR(rms[0]);
        free(ionos);free(rates);free(rms);free(quals);
        free(ionos_);free(rates_);free(sqrms);free(sqrms_);free(quals_);
        return 1;
    } else {
        for (i=0;i<n;i++) { /* for each grid */
            /* search stec data */
            if (!stec_data(nav->stec+index[i],time,sat,ionos+i,rates+i,rms+i,quals+i,&slip)) {
                free(ionos);free(rates);free(rms);free(quals);
                free(ionos_);free(rates_);free(sqrms);free(sqrms_);free(quals_);
                return 0;
            }
            if (slip) *brk=1;
        }

        *iono=*rate=*var=0.0;

        if (n==4&&Gmat&&Emat) {
            /* Emat, Gmat */
            for (i=0;i<n;i++) sqrms[i]=SQR(rms[i]);
            matmul("NN",n,1,n,1.0,Gmat,ionos,0.0,ionos_);
            matmul("NN",n,1,n,1.0,Gmat,rates,0.0,rates_);
            matmul("NN",n,1,n,1.0,Gmat,sqrms,0.0,sqrms_);
            matmul("NN",n,1,n,1.0,Gmat,quals,0.0,quals_);
            *iono=dot(Emat,ionos_,4);
            *rate=dot(Emat,rates_,4);
            *var =dot(Emat,sqrms_,4);
            *quality=dot(Emat,quals_,4);
        } else {
            /* weight */
            for (i=0;i<n;i++) {
                *iono+=ionos[i]*weight[i];
                *rate+=rates[i]*weight[i];
                *var+=SQR(rms[i])*weight[i];
                *quality+=quals[i]*weight[i];
            }
        }
    }

    free(ionos);free(rates);free(rms);free(quals);
    free(ionos_);free(rates_);free(sqrms);free(sqrms_);free(quals_);
    
    return 1;
}

extern int trop_grid_data(const nav_t *nav, const int *index, gtime_t time,
                          int n, const double *weight, const double *Gmat,
                          const double *Emat, double *zwd, double *ztd,
                          double *quality, int *tbrk)
{
    int i,*valid,valid_=0;
    double *zwds,*ztds,*zwds_,*ztds_,*quals,*quals_;
    
    if (n<=0) return 0;

    zwds=mat(n,1);ztds=mat(n,1);quals=mat(n,1);valid=imat(n,1);

    if (n == 1) {
        if (!trop_data(nav->zwd+index[0],time,ztd,zwd,quals,index[0],valid)) {
            free(zwds);free(ztds);free(quals);free(valid);
            return 0;
        }
        return 1;
    } else {
        for (i=0;i<n;i++) { /* for each grid */
            /* search stec data */
            if (!trop_data(nav->zwd+index[i],time,ztds+i,zwds+i,quals+i,index[i],valid+i)) {
                free(zwds);free(ztds);free(quals);free(valid);
                return 0;
            }
        }
    }

    *zwd=*ztd=*quality=0.0;

    if (n==4&&Gmat&&Emat) {
        zwds_=mat(n,1);ztds_=mat(n,1);quals_=mat(n,1);
        /* Emat, Gmat */
        matmul("NN",n,1,n,1.0,Gmat,ztds,0.0,ztds_);
        matmul("NN",n,1,n,1.0,Gmat,zwds,0.0,zwds_);
        matmul("NN",n,1,n,1.0,Gmat,quals,0.0,quals_);
        *ztd=dot(Emat,ztds_,4);
        *zwd=dot(Emat,zwds_,4);
        *quality=dot(Emat,quals_,4);
        free(zwds_);free(ztds_);free(quals_);
    } else {
        /* weight */
        for (i=0;i<n;i++) {
            *zwd+=zwds[i]*weight[i];
            *ztd+=ztds[i]*weight[i];
            *quality+=quals[i]*weight[i];
        }
    }
    for (i=0;i<n;i++) valid_+=valid[i];
    *tbrk=valid_==n?0:1;
    
    free(zwds);free(ztds);free(quals);free(valid);
    
    return 1;
}
extern void get_oload_ptr(nav_t *nav) {
    int i;
    for (i=0; i<CSSR_MAX_NETWORK; i++) {
        nav->oload[i] =(olod_t *)&clas_oload[i];
        nav->oload[i]->gridnum=0;
    }
}