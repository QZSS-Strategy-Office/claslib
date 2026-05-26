/*------------------------------------------------------------------------------
* cssr.h : Compact SSR constants, types and function prototypes
*
*          Copyright (C) 2015- by Mitsubishi Electric Corporation, All rights reserved.
*-----------------------------------------------------------------------------*/
#ifndef CSSR_H
#define CSSR_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#if !defined(__CYGWIN__)
#include <stdint.h>
#endif

/* constants -----------------------------------------------------------------*/

#define VER_CSSR  "0.2"             /* library version */

/* cssr */
#define CSSR_MAX_GNSS     16
#define CSSR_MAX_SV_GNSS  40
#define CSSR_MAX_SV       64
#define CSSR_MAX_SIG      16
#define CSSR_MAX_CELLMASK 64
#define CSSR_MAX_NET      32
#define CSSR_MAX_LOCAL_SV 32
#define CSSR_MAX_GP       128
#define CSSR_MAX_NETWORK  32

#define CSSR_SYS_GPS    0
#define CSSR_SYS_GLO    1
#define CSSR_SYS_GAL    2
#define CSSR_SYS_BDS    3
#define CSSR_SYS_QZS    4
#define CSSR_SYS_SBS    5
#define CSSR_SYS_NONE   -1

#define CSSR_TYPE_NUM   14
#define CSSR_TYPE_MASK  1
#define CSSR_TYPE_OC    2
#define CSSR_TYPE_CC    3
#define CSSR_TYPE_CB    4
#define CSSR_TYPE_PB    5
#define CSSR_TYPE_BIAS  6
#define CSSR_TYPE_URA   7
#define CSSR_TYPE_STEC  8
#define CSSR_TYPE_GRID  9
#define CSSR_TYPE_SI    10
#define CSSR_TYPE_COMBO 11
#define CSSR_TYPE_ATMOS 12

#define CSSR_TYPE_INIT  254
#define CSSR_TYPE_NULL  255

#define P2_S16_MAX 32767
#define P2_S15_MAX 16383
#define P2_S14_MAX 8191
#define P2_S13_MAX 4095
#define P2_S12_MAX 2047
#define P2_S11_MAX 1023
#define P2_S10_MAX 511
#define P2_S9_MAX  255
#define P2_S8_MAX  127
#define P2_S7_MAX  63
#define P2_S6_MAX  31

#define CSSR_TROP_HS_REF    2.3
#define CSSR_TROP_WET_REF   0.252

#define CSSR_UPDATE_TROP	0
#define CSSR_UPDATE_STEC	1
#define CSSR_UPDATE_PBIAS	2

#define TROPVALIDAGE  3600
#define STECVALIDAGE  3600

#define INVALID_VALUE -10000

typedef struct {
    int stec_type;
} cssropt_t;

typedef struct {
    gtime_t t0[2];
    double udi[2];
    int iod[2];
    int ngp;
    float quality_f[CSSR_MAX_LOCAL_SV];
    float trop_wet[CSSR_MAX_GP];
    float trop_total[CSSR_MAX_GP];
    int nsat_f;
    int sat_f[CSSR_MAX_LOCAL_SV];
    float quality;
    double a[CSSR_MAX_LOCAL_SV][4];
    int nsat[CSSR_MAX_GP];
    int sat[CSSR_MAX_GP][CSSR_MAX_LOCAL_SV];
    float stec[CSSR_MAX_GP][CSSR_MAX_LOCAL_SV];
    double grid[CSSR_MAX_GP][3];
    int update[3];
} ssrn_t;

typedef struct {
    int ver;
    cssropt_t opt;
    int iod;
    int iod_sv;
    int inet;
    int week;
    uint8_t cmi[CSSR_MAX_GNSS]; /* cellmask existence flag */
    uint64_t svmask[CSSR_MAX_GNSS];
    uint16_t sigmask[CSSR_MAX_GNSS];
    uint16_t cellmask[CSSR_MAX_SV];
    uint64_t net_svmask[CSSR_MAX_NET];
    int ngnss;
    int nsat;
    int ncell;
    int sat[CSSR_MAX_SV];
    int nsat_n[CSSR_MAX_NET];
    int sat_n[CSSR_MAX_NET][CSSR_MAX_LOCAL_SV];
    int nsig[CSSR_MAX_SV];
    int sigmask_s[CSSR_MAX_SV];
    int amb_bias[MAXSAT][MAXCODE];
    uint8_t disc[MAXSAT][MAXCODE];
    float quality_i;    /* ionosphere quality */
    int l6delivery;
    int l6facility;
    ssrn_t ssrn[CSSR_MAX_NET];
    int si_cnt;
    int si_sz;
    uint64_t si_data[4];
} cssr_t;

extern FILE *open_L6(char **infile, int n);
extern int dumpcssr(char **infile, int n, FILE **ofp, filopt_t *fopt);
extern int open_outputfiles(FILE **ofp);
extern void close_outputfiles(FILE **ofp);

#endif /* CSSR_H */
