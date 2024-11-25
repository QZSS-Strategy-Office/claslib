/*------------------------------------------------------------------------------
* cssr.h : Compact SSR Grid Definition
*
*          Copyright (C) 2015- by Mitsubishi Electric Corporation, All rights reserved.
*-----------------------------------------------------------------------------*/
#ifndef SRC_CLASGRID_H_
#define SRC_CLASGRID_H_

#ifndef CLASGRID_GLOBAL_DEFINE

extern double clas_grid[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][3];
extern int grid_stat[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
extern olod_t clas_oload[CSSR_MAX_NETWORK];

#else

int grid_stat[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP];
olod_t clas_oload[CSSR_MAX_NETWORK];
double clas_grid[CSSR_MAX_NETWORK][RTCM_SSR_MAX_GP][3];

#endif /* CLASGRID_GLOBAL_DEFINE */

#endif /* SRC_CLASGRID_H_ */
