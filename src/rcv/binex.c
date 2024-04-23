/*------------------------------------------------------------------------------
* binex.c : binex dependent functions
*
*          Copyright (C) 2013-2018 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] UNAVCO, BINEX: Binary exchange format
*         (http://binex.unavco.org/binex.html)
*
* version : $Revision:$ $Date:$
* history : 2013/02/20 1.0 new
*           2013/04/15 1.1 support 0x01-05 beidou-2/compass ephemeris
*           2013/05/18 1.2 fix bug on decoding obsflags in message 0x7f-05
*           2014/04/27 1.3 fix bug on decoding iode for message 0x01-02
*           2015/12/05 1.4 fix bug on decoding tgd for message 0x01-05
*           2016/07/29 1.5 crc16() -> rtk_crc16()
*           2017/04/11 1.6 (char *) -> (signed char *)
*                          fix bug on unchange-test of beidou ephemeris
*           2017/11/28 1.7 fix bug on decoding galileo ephemeris (0x01-04)
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include <string.h>

#define BNXSYNC1    0xC2    /* binex sync (little-endian,regular-crc) */
#define BNXSYNC2    0xE2    /* binex sync (big-endian   ,regular-crc) */
#define BNXSYNC3    0xC8    /* binex sync (little-endian,enhanced-crc) */
#define BNXSYNC4    0xE8    /* binex sync (big-endian   ,enhanced-crc) */

#define BNXSYNC1R   0xD2    /* binex sync (little-endian,regular-crc,rev) */
#define BNXSYNC2R   0xF2    /* binex sync (big-endian   ,regular-crc,rev) */
#define BNXSYNC3R   0xD8    /* binex sync (little-endian,enhanced-crc,rev) */
#define BNXSYNC4R   0xF8    /* binex sync (big-endian   ,enhanced-crc,rev) */

#define MIN(x,y)    ((x)<(y)?(x):(y))

/* ura table -----------------------------------------------------------------*/
static const double ura_eph[]={
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
    3072.0,6144.0,0.0
};
/* get fields (big-endian) ---------------------------------------------------*/
#define U1(p) (*((unsigned char *)(p)))
#define I1(p) (*((signed char *)(p)))

static unsigned short U2(unsigned char *p)
{
    unsigned short value;
    unsigned char *q=(unsigned char *)&value+1;
    int i;
    for (i=0;i<2;i++) *q--=*p++;
    return value;
}
static unsigned int U4(unsigned char *p)
{
    unsigned int value;
    unsigned char *q=(unsigned char *)&value+3;
    int i;
    for (i=0;i<4;i++) *q--=*p++;
    return value;
}
static int I4(unsigned char *p)
{
    return (int)U4(p);
}
static float R4(unsigned char *p)
{
    float value;
    unsigned char *q=(unsigned char *)&value+3;
    int i;
    for (i=0;i<4;i++) *q--=*p++;
    return value;
}
static double R8(unsigned char *p)
{
    double value;
    unsigned char *q=(unsigned char *)&value+7;
    int i;
    for (i=0;i<8;i++) *q--=*p++;
    return value;
}
/* get binex 1-4 byte unsigned integer (big endian) --------------------------*/
static int getbnxi(unsigned char *p, unsigned int *val)
{
    int i;
    
    for (*val=0,i=0;i<3;i++) {
        *val=(*val<<7)+(p[i]&0x7F);
        if (!(p[i]&0x80)) return i+1;
    }
    *val=(*val<<8)+p[i];
    return 4;
}
/* checksum 8 parity ---------------------------------------------------------*/
static unsigned char csum8(const unsigned char *buff, int len)
{
    unsigned char cs=0;
    int i;
    
    for (i=0;i<len;i++) {
        cs^=buff[i];
    }
    return cs;
}
/* adjust weekly rollover of gps time ----------------------------------------*/
static gtime_t adjweek(gtime_t time, double tow)
{
    double tow_p;
    int week;
    tow_p=time2gpst(time,&week);
    if      (tow<tow_p-302400.0) tow+=604800.0;
    else if (tow>tow_p+302400.0) tow-=604800.0;
    return gpst2time(week,tow);
}
/* adjust daily rollover of time ---------------------------------------------*/
static gtime_t adjday(gtime_t time, double tod)
{
    double ep[6],tod_p;
    time2epoch(time,ep);
    tod_p=ep[3]*3600.0+ep[4]*60.0+ep[5];
    if      (tod<tod_p-43200.0) tod+=86400.0;
    else if (tod>tod_p+43200.0) tod-=86400.0;
    ep[3]=ep[4]=ep[5]=0.0;
    return timeadd(epoch2time(ep),tod);
}
/* ura value (m) to ura index ------------------------------------------------*/
static int uraindex(double value)
{
    int i;
    for (i=0;i<15;i++) if (ura_eph[i]>=value) break;
    return i;
}
/* decode binex mesaage 0x00-00: comment -------------------------------------*/
static int decode_bnx_00_00(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-00: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-01: program or software package -----------------*/
static int decode_bnx_00_01(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-01: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-02: program operator ----------------------------*/
static int decode_bnx_00_02(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-02: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-03: reserved ------------------------------------*/
static int decode_bnx_00_03(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-03: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-04: site name/description -----------------------*/
static int decode_bnx_00_04(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-04: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-05: site number ---------------------------------*/
static int decode_bnx_00_05(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-05: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-06: monumnent name ------------------------------*/
static int decode_bnx_00_06(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-06: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-07: monumnent number ----------------------------*/
static int decode_bnx_00_07(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-07: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-08: marker name ---------------------------------*/
static int decode_bnx_00_08(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-08: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-09: marker number -------------------------------*/
static int decode_bnx_00_09(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-09: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-0a: reference point name ------------------------*/
static int decode_bnx_00_0a(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-0a: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-0b: reference point number ----------------------*/
static int decode_bnx_00_0b(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-0b: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-0c: date esttablished ---------------------------*/
static int decode_bnx_00_0c(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-0c: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-0d: reserved ------------------------------------*/
static int decode_bnx_00_0d(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-0d: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-0e: reserved ------------------------------------*/
static int decode_bnx_00_0e(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-0e: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-0f: 4-character id ------------------------------*/
static int decode_bnx_00_0f(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-0f: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-10: project name --------------------------------*/
static int decode_bnx_00_10(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-10: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-11: principal investigator for this project -----*/
static int decode_bnx_00_11(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-11: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-12: pi's agency/institution ---------------------*/
static int decode_bnx_00_12(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-12: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-13: pi's contact information --------------------*/
static int decode_bnx_00_13(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-13: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-14: site operator -------------------------------*/
static int decode_bnx_00_14(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-14: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-15: site operator's agency/institution ----------*/
static int decode_bnx_00_15(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-15: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-16: site operator's contact information ---------*/
static int decode_bnx_00_16(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-16: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-17: antenna type --------------------------------*/
static int decode_bnx_00_17(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-17: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-18: antenna number ------------------------------*/
static int decode_bnx_00_18(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-18: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-19: receiver type -------------------------------*/
static int decode_bnx_00_19(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-19: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-1a: receiver number -----------------------------*/
static int decode_bnx_00_1a(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-1a: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-1b: receiver firmware version -------------------*/
static int decode_bnx_00_1b(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-1b: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-1c: antenna mount description -------------------*/
static int decode_bnx_00_1c(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-1c: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-1d: antenna xyz position ------------------------*/
static int decode_bnx_00_1d(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-1d: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-1e: antenna geographic position -----------------*/
static int decode_bnx_00_1e(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-1e: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-1f: antenna offset from reference point ---------*/
static int decode_bnx_00_1f(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-1f: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-20: antenna radome type -------------------------*/
static int decode_bnx_00_20(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-20: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-21: antenna radome number -----------------------*/
static int decode_bnx_00_21(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-21: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-22: geocode -------------------------------------*/
static int decode_bnx_00_22(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-22: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00-7f: notes/additional information ----------------*/
static int decode_bnx_00_7f(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x00-7f: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x00: site/monument/marker/ref point/setup metadata --*/
static int decode_bnx_00(raw_t *raw, unsigned char *buff, int len)
{
    const static double gpst0[]={1980,1,6,0,0,0};
    char *msg;
    unsigned char *p=buff;
    unsigned int min,qsec,src,fid;
    int n=6;
    
    min =U4(p); p+=4;
    qsec=U1(p); p+=1;
    src =U1(p); p+=1;
    n+=getbnxi(p,&fid);
    raw->time=timeadd(epoch2time(gpst0),min*60.0+qsec*0.25);
    
    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg," fid=%02X time=%s src=%d",fid,time_str(raw->time,0),src);
    }
    switch (fid) {
        case 0x00: return decode_bnx_00_00(raw,buff+n,len-n);
        case 0x01: return decode_bnx_00_01(raw,buff+n,len-n);
        case 0x02: return decode_bnx_00_02(raw,buff+n,len-n);
        case 0x03: return decode_bnx_00_03(raw,buff+n,len-n);
        case 0x04: return decode_bnx_00_04(raw,buff+n,len-n);
        case 0x05: return decode_bnx_00_05(raw,buff+n,len-n);
        case 0x06: return decode_bnx_00_06(raw,buff+n,len-n);
        case 0x07: return decode_bnx_00_07(raw,buff+n,len-n);
        case 0x08: return decode_bnx_00_08(raw,buff+n,len-n);
        case 0x09: return decode_bnx_00_09(raw,buff+n,len-n);
        case 0x0A: return decode_bnx_00_0a(raw,buff+n,len-n);
        case 0x0B: return decode_bnx_00_0b(raw,buff+n,len-n);
        case 0x0C: return decode_bnx_00_0c(raw,buff+n,len-n);
        case 0x0D: return decode_bnx_00_0d(raw,buff+n,len-n);
        case 0x0E: return decode_bnx_00_0e(raw,buff+n,len-n);
        case 0x0F: return decode_bnx_00_0f(raw,buff+n,len-n);
        case 0x10: return decode_bnx_00_10(raw,buff+n,len-n);
        case 0x11: return decode_bnx_00_11(raw,buff+n,len-n);
        case 0x12: return decode_bnx_00_12(raw,buff+n,len-n);
        case 0x13: return decode_bnx_00_13(raw,buff+n,len-n);
        case 0x14: return decode_bnx_00_14(raw,buff+n,len-n);
        case 0x15: return decode_bnx_00_15(raw,buff+n,len-n);
        case 0x16: return decode_bnx_00_16(raw,buff+n,len-n);
        case 0x17: return decode_bnx_00_17(raw,buff+n,len-n);
        case 0x18: return decode_bnx_00_18(raw,buff+n,len-n);
        case 0x19: return decode_bnx_00_19(raw,buff+n,len-n);
        case 0x1A: return decode_bnx_00_1a(raw,buff+n,len-n);
        case 0x1B: return decode_bnx_00_1b(raw,buff+n,len-n);
        case 0x1C: return decode_bnx_00_1c(raw,buff+n,len-n);
        case 0x1D: return decode_bnx_00_1d(raw,buff+n,len-n);
        case 0x1E: return decode_bnx_00_1e(raw,buff+n,len-n);
        case 0x1F: return decode_bnx_00_1f(raw,buff+n,len-n);
        case 0x20: return decode_bnx_00_20(raw,buff+n,len-n);
        case 0x21: return decode_bnx_00_21(raw,buff+n,len-n);
        case 0x22: return decode_bnx_00_22(raw,buff+n,len-n);
        case 0x7F: return decode_bnx_00_7f(raw,buff+n,len-n);
    }
    return 0;
}
/* decode binex mesaage 0x01-00: coded (raw bytes) gnss ephemeris ------------*/
static int decode_bnx_01_00(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x01-00: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x01-01: decoded gps ephmemeris ----------------------*/
static int decode_bnx_01_01(raw_t *raw, unsigned char *buff, int len)
{
    eph_t eph={0};
    unsigned char *p=buff;
    double tow,ura,sqrtA;
    int prn,flag;
    
    trace(4,"binex 0x01-01: len=%d\n",len);
    
    if (len>=127) {
        prn       =U1(p)+1;      p+=1;
        eph.week  =U2(p);        p+=2;
        tow       =I4(p);        p+=4;
        eph.toes  =I4(p);        p+=4;
        eph.tgd[0]=R4(p);        p+=4;
        eph.iodc  =I4(p);        p+=4;
        eph.f2    =R4(p);        p+=4;
        eph.f1    =R4(p);        p+=4;
        eph.f0    =R4(p);        p+=4;
        eph.iode  =I4(p);        p+=4;
        eph.deln  =R4(p)*SC2RAD; p+=4;
        eph.M0    =R8(p);        p+=8;
        eph.e     =R8(p);        p+=8;
        sqrtA     =R8(p);        p+=8;
        eph.cic   =R4(p);        p+=4;
        eph.crc   =R4(p);        p+=4;
        eph.cis   =R4(p);        p+=4;
        eph.crs   =R4(p);        p+=4;
        eph.cuc   =R4(p);        p+=4;
        eph.cus   =R4(p);        p+=4;
        eph.OMG0  =R8(p);        p+=8;
        eph.omg   =R8(p);        p+=8;
        eph.i0    =R8(p);        p+=8;
        eph.OMGd  =R4(p)*SC2RAD; p+=4;
        eph.idot  =R4(p)*SC2RAD; p+=4;
        ura       =R4(p)*0.1;    p+=4;
        eph.svh   =U2(p);        p+=2;
        flag      =U2(p);
    }
    else {
        trace(2,"binex 0x01-01: length error len=%d\n",len);
        return -1;
    }
    if (!(eph.sat=satno(SYS_GPS,prn))) {
        trace(2,"binex 0x01-01: satellite error prn=%d\n",prn);
        return -1;
    }
    eph.A=sqrtA*sqrtA;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,eph.toes);
    eph.ttr=adjweek(eph.toe,tow);
    eph.fit=flag&0xFF;
    eph.flag=(flag>>8)&0x01;
    eph.code=(flag>>9)&0x03;
    eph.sva=ura_index(ura,SYS_GPS);
    
    if (!strstr(raw->opt,"-EPHALL")) {
        if (raw->nav.eph[eph.sat-1].iode==eph.iode&&
            raw->nav.eph[eph.sat-1].iodc==eph.iodc) return 0; /* unchanged */
    }
    raw->nav.eph[eph.sat-1]=eph;
    raw->ephsat=eph.sat;
    return 2;
}
/* decode binex mesaage 0x01-02: decoded glonass ephmemeris ------------------*/
static int decode_bnx_01_02(raw_t *raw, unsigned char *buff, int len)
{
    geph_t geph={0};
    unsigned char *p=buff;
    double tod,tof,tau_gps;
    int prn,day,leap;
    
    trace(4,"binex 0x01-02: len=%d\n",len);
    
    if (len>=119) {
        prn        =U1(p)+1;   p+=1;
        day        =U2(p);     p+=2;
        tod        =U4(p);     p+=4;
        geph.taun  =-R8(p);    p+=8;
        geph.gamn  =R8(p);     p+=8;
        tof        =U4(p);     p+=4;
        geph.pos[0]=R8(p)*1E3; p+=8;
        geph.vel[0]=R8(p)*1E3; p+=8;
        geph.acc[0]=R8(p)*1E3; p+=8;
        geph.pos[1]=R8(p)*1E3; p+=8;
        geph.vel[1]=R8(p)*1E3; p+=8;
        geph.acc[1]=R8(p)*1E3; p+=8;
        geph.pos[2]=R8(p)*1E3; p+=8;
        geph.vel[2]=R8(p)*1E3; p+=8;
        geph.acc[2]=R8(p)*1E3; p+=8;
        geph.svh   =U1(p)&0x1; p+=1;
        geph.frq   =I1(p);     p+=1;
        geph.age   =U1(p);     p+=1;
        leap       =U1(p);     p+=1;
        tau_gps    =R8(p);     p+=8;
        geph.dtaun =R8(p);
    }
    else {
        trace(2,"binex 0x01-02: length error len=%d\n",len);
        return -1;
    }
    if (!(geph.sat=satno(SYS_GLO,prn))) {
        trace(2,"binex 0x01-02: satellite error prn=%d\n",prn);
        return -1;
    }
    if (raw->time.time==0) return 0;
    geph.toe=utc2gpst(adjday(raw->time,tod-10800.0));
    geph.tof=utc2gpst(adjday(raw->time,tof-10800.0));
    geph.iode=(int)(fmod(tod,86400.0)/900.0+0.5);
    
    if (!strstr(raw->opt,"-EPHALL")) {
        if (fabs(timediff(geph.toe,raw->nav.geph[prn-MINPRNGLO].toe))<1.0&&
            geph.svh==raw->nav.geph[prn-MINPRNGLO].svh) return 0; /* unchanged */
    }
    raw->nav.geph[prn-1]=geph;
    raw->ephsat=geph.sat;
    return 2;
}
/* decode binex mesaage 0x01-03: decoded sbas ephmemeris ---------------------*/
static int decode_bnx_01_03(raw_t *raw, unsigned char *buff, int len)
{
    seph_t seph={0};
    unsigned char *p=buff;
    double tow,tod,tof;
    int prn,week,iodn;
    
    trace(4,"binex 0x01-03: len=%d\n",len);
    
    if (len>=98) {
        prn        =U1(p);     p+=1;
        week       =U2(p);     p+=2;
        tow        =U4(p);     p+=4;
        seph.af0   =R8(p);     p+=8;
        tod        =R4(p);     p+=4;
        tof        =U4(p);     p+=4;
        seph.pos[0]=R8(p)*1E3; p+=8;
        seph.vel[0]=R8(p)*1E3; p+=8;
        seph.acc[0]=R8(p)*1E3; p+=8;
        seph.pos[1]=R8(p)*1E3; p+=8;
        seph.vel[1]=R8(p)*1E3; p+=8;
        seph.acc[1]=R8(p)*1E3; p+=8;
        seph.pos[2]=R8(p)*1E3; p+=8;
        seph.vel[2]=R8(p)*1E3; p+=8;
        seph.acc[2]=R8(p)*1E3; p+=8;
        seph.svh   =U1(p);     p+=1;
        seph.sva   =U1(p);     p+=1;
        iodn       =U1(p);
    }
    else {
        trace(2,"binex 0x01-03 length error: len=%d\n",len);
        return -1;
    }
    if (!(seph.sat=satno(SYS_SBS,prn))) {
        trace(2,"binex 0x01-03 satellite error: prn=%d\n",prn);
        return -1;
    }
    seph.t0=gpst2time(week,tow);
    seph.tof=adjweek(seph.t0,tof);
    
    if (!strstr(raw->opt,"-EPHALL")) {
        if (fabs(timediff(seph.t0,raw->nav.seph[prn-MINPRNSBS].t0))<1.0&&
            seph.sva==raw->nav.seph[prn-MINPRNSBS].sva) return 0; /* unchanged */
    }
    raw->nav.seph[prn-MINPRNSBS]=seph;
    raw->ephsat=seph.sat;
    return 2;
}
/* decode binex mesaage 0x01-04: decoded galileo ephmemeris ------------------*/
static int decode_bnx_01_04(raw_t *raw, unsigned char *buff, int len)
{
    eph_t eph={0};
    unsigned char *p=buff;
    double tow,ura,sqrtA;
    int prn;
    
    trace(4,"binex 0x01-04: len=%d\n",len);
    
    if (len>=127) {
        prn       =U1(p)+1;      p+=1;
        eph.week  =U2(p);        p+=2; /* gal-week = gps-week */
        tow       =I4(p);        p+=4;
        eph.toes  =I4(p);        p+=4;
        eph.tgd[0]=R4(p);        p+=4; /* BGD E5a/E1 */
        eph.tgd[1]=R4(p);        p+=4; /* BGD E5b/E1 */
        eph.iode  =I4(p);        p+=4; /* IODnav */
        eph.f2    =R4(p);        p+=4;
        eph.f1    =R4(p);        p+=4;
        eph.f0    =R4(p);        p+=4;
        eph.deln  =R4(p)*SC2RAD; p+=4;
        eph.M0    =R8(p);        p+=8;
        eph.e     =R8(p);        p+=8;
        sqrtA     =R8(p);        p+=8;
        eph.cic   =R4(p);        p+=4;
        eph.crc   =R4(p);        p+=4;
        eph.cis   =R4(p);        p+=4;
        eph.crs   =R4(p);        p+=4;
        eph.cuc   =R4(p);        p+=4;
        eph.cus   =R4(p);        p+=4;
        eph.OMG0  =R8(p);        p+=8;
        eph.omg   =R8(p);        p+=8;
        eph.i0    =R8(p);        p+=8;
        eph.OMGd  =R4(p)*SC2RAD; p+=4;
        eph.idot  =R4(p)*SC2RAD; p+=4;
        ura       =R4(p);        p+=4;
        eph.svh   =U2(p);        p+=2;
        eph.code  =U2(p);              /* data source */
    }
    else {
        trace(2,"binex 0x01-04: length error len=%d\n",len);
        return -1;
    }
    if (!(eph.sat=satno(SYS_GAL,prn))) {
        trace(2,"binex 0x01-04: satellite error prn=%d\n",prn);
        return -1;
    }
    eph.A=sqrtA*sqrtA;
    eph.iodc=eph.iode;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,eph.toes);
    eph.ttr=adjweek(eph.toe,tow);
    if (ura < 0) {
        /*  -(SISA index + 1) for SISA index = 0-255 */
        eph.sva=-(ura + 1);
    }
    else {
        /* accuracy in distance in meters for SISA index = 0-125 */
        eph.sva=ura_index(ura,SYS_GAL);
    }

    if (!strstr(raw->opt,"-EPHALL")) {
        if (raw->nav.eph[eph.sat-1].iode==eph.iode&&
            raw->nav.eph[eph.sat-1].iodc==eph.iodc) return 0; /* unchanged */
    }
    raw->nav.eph[eph.sat-1]=eph;
    raw->ephsat=eph.sat;
    return 2;
}
/* beidou signed 10 bit tgd -> sec -------------------------------------------*/
static double bds_tgd(int tgd)
{
    tgd&=0x3FF;
    return (tgd&0x200)?-1E-10*((~tgd)&0x1FF):1E-10*(tgd&0x1FF);
}
/* decode binex mesaage 0x01-05: decoded beidou-2/compass ephmemeris ---------*/
static int decode_bnx_01_05(raw_t *raw, unsigned char *buff, int len)
{
    eph_t eph={0};
    unsigned char *p=buff;
    double tow,toc,sqrtA;
    int prn,flag1,flag2;
    
    trace(4,"binex 0x01-05: len=%d\n",len);
    
    if (len>=117) {
        prn       =U1(p);        p+=1;
        eph.week  =U2(p);        p+=2;
        tow       =I4(p);        p+=4;
        toc       =I4(p);        p+=4;
        eph.toes  =I4(p);        p+=4;
        eph.f2    =R4(p);        p+=4;
        eph.f1    =R4(p);        p+=4;
        eph.f0    =R4(p);        p+=4;
        eph.deln  =R4(p)*SC2RAD; p+=4;
        eph.M0    =R8(p);        p+=8;
        eph.e     =R8(p);        p+=8;
        sqrtA     =R8(p);        p+=8;
        eph.cic   =R4(p);        p+=4;
        eph.crc   =R4(p);        p+=4;
        eph.cis   =R4(p);        p+=4;
        eph.crs   =R4(p);        p+=4;
        eph.cuc   =R4(p);        p+=4;
        eph.cus   =R4(p);        p+=4;
        eph.OMG0  =R8(p);        p+=8;
        eph.omg   =R8(p);        p+=8;
        eph.i0    =R8(p);        p+=8;
        eph.OMGd  =R4(p)*SC2RAD; p+=4;
        eph.idot  =R4(p)*SC2RAD; p+=4;
        flag1     =U2(p);        p+=2;
        flag2     =U4(p);
    }
    else {
        trace(2,"binex 0x01-05: length error len=%d\n",len);
        return -1;
    }
    if (!(eph.sat=satno(SYS_CMP,prn))) {
        trace(2,"binex 0x01-05: satellite error prn=%d\n",prn);
        return 0;
    }
    eph.A=sqrtA*sqrtA;
    eph.toe=gpst2time(eph.week+1356,eph.toes+14.0); /* bdt -> gpst */
    eph.toc=gpst2time(eph.week+1356,eph.toes+14.0); /* bdt -> gpst */
    eph.ttr=adjweek(eph.toe,tow+14.0); /* bdt -> gpst */
    eph.iodc=(flag1>>1)&0x1F;
    eph.iode=(flag1>>6)&0x1F;
    eph.svh=flag1&0x01;
    eph.sva=flag2&0x0F; /* ura index */
    eph.tgd[0]=bds_tgd(flag2>> 4); /* TGD1 (s) */
    eph.tgd[1]=bds_tgd(flag2>>14); /* TGD2 (s) */
    eph.flag=(flag1>>11)&0x07; /* nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    eph.code=(flag2>>25)&0x7F;
        /* message source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q)*/
    
    if (!strstr(raw->opt,"-EPHALL")) {
        if (timediff(raw->nav.eph[eph.sat-1].toe,eph.toe)==0.0&&
            raw->nav.eph[eph.sat-1].iode==eph.iode&&
            raw->nav.eph[eph.sat-1].iodc==eph.iodc) return 0; /* unchanged */
    }
    raw->nav.eph[eph.sat-1]=eph;
    raw->ephsat=eph.sat;
    return 2;
}
/* decode binex mesaage 0x01-06: decoded qzss ephmemeris ---------------------*/
static int decode_bnx_01_06(raw_t *raw, unsigned char *buff, int len)
{
    eph_t eph={0};
    unsigned char *p=buff;
    double tow,ura,sqrtA;
    int prn,flag;
    
    trace(4,"binex 0x01-06: len=%d\n",len);
    
    if (len>=127) {
        prn       =U1(p);        p+=1;
        eph.week  =U2(p);        p+=2;
        tow       =I4(p);        p+=4;
        eph.toes  =I4(p);        p+=4;
        eph.tgd[0]=R4(p);        p+=4;
        eph.iodc  =I4(p);        p+=4;
        eph.f2    =R4(p);        p+=4;
        eph.f1    =R4(p);        p+=4;
        eph.f0    =R4(p);        p+=4;
        eph.iode  =I4(p);        p+=4;
        eph.deln  =R4(p)*SC2RAD; p+=4;
        eph.M0    =R8(p);        p+=8;
        eph.e     =R8(p);        p+=8;
        sqrtA     =R8(p);        p+=8;
        eph.cic   =R4(p);        p+=4;
        eph.crc   =R4(p);        p+=4;
        eph.cis   =R4(p);        p+=4;
        eph.crs   =R4(p);        p+=4;
        eph.cuc   =R4(p);        p+=4;
        eph.cus   =R4(p);        p+=4;
        eph.OMG0  =R8(p);        p+=8;
        eph.omg   =R8(p);        p+=8;
        eph.i0    =R8(p);        p+=8;
        eph.OMGd  =R4(p)*SC2RAD; p+=4;
        eph.idot  =R4(p)*SC2RAD; p+=4;
        ura       =R4(p)*0.1;    p+=4;
        eph.svh   =U2(p);        p+=2;
        flag      =U2(p);
    }
    else {
        trace(2,"binex 0x01-06: length error len=%d\n",len);
        return -1;
    }
    if (!(eph.sat=satno(SYS_QZS,prn))) {
        trace(2,"binex 0x01-06: satellite error prn=%d\n",prn);
        return 0;
    }
    eph.A=sqrtA*sqrtA;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,eph.toes);
    eph.ttr=adjweek(eph.toe,tow);
    eph.fit=(flag&0x01)?0.0:2.0; /* 0:2hr,1:>2hr */
    eph.sva=ura_index(ura,SYS_QZS);
    eph.code=2; /* codes on L2 channel */
    
    if (!strstr(raw->opt,"-EPHALL")) {
        if (raw->nav.eph[eph.sat-1].iode==eph.iode&&
            raw->nav.eph[eph.sat-1].iodc==eph.iodc) return 0; /* unchanged */
    }
    raw->nav.eph[eph.sat-1]=eph;
    raw->ephsat=eph.sat;
    return 2;
}
/* decode binex mesaage 0x01-14: upgraded decoded galileo ephemeris ----------*/
static int decode_bnx_01_14(raw_t *raw, unsigned char *buff, int len)
{
    eph_t eph={0};
    unsigned char *p=buff;
    double tow,ura,sqrtA;
    int prn,tocs;
    
    trace(4,"binex 0x01-14: len=%d\n",len);
    
    if (len>=135) {
        prn       =U1(p)+1;      p+=1;
        eph.week  =U2(p);        p+=2; /* gal-week = gps-week */
        tow       =I4(p);        p+=4;
        tocs      =I4(p);        p+=4;
        eph.toes  =I4(p);        p+=4;
        eph.tgd[0]=R4(p);        p+=4; /* BGD E5a/E1 */
        eph.tgd[1]=R4(p);        p+=4; /* BGD E5b/E1 */
        eph.iode  =I4(p);        p+=4; /* IODnav */
        eph.f2    =R4(p);        p+=4;
        eph.f1    =R4(p);        p+=4;
        eph.f0    =R8(p);        p+=8;
        eph.deln  =R4(p)*SC2RAD; p+=4;
        eph.M0    =R8(p);        p+=8;
        eph.e     =R8(p);        p+=8;
        sqrtA     =R8(p);        p+=8;
        eph.cic   =R4(p);        p+=4;
        eph.crc   =R4(p);        p+=4;
        eph.cis   =R4(p);        p+=4;
        eph.crs   =R4(p);        p+=4;
        eph.cuc   =R4(p);        p+=4;
        eph.cus   =R4(p);        p+=4;
        eph.OMG0  =R8(p);        p+=8;
        eph.omg   =R8(p);        p+=8;
        eph.i0    =R8(p);        p+=8;
        eph.OMGd  =R4(p)*SC2RAD; p+=4;
        eph.idot  =R4(p)*SC2RAD; p+=4;
        ura       =R4(p);        p+=4;
        eph.svh   =U2(p);        p+=2;
        eph.code  =U2(p);              /* data source */
    }
    else {
        trace(2,"binex 0x01-14: length error len=%d\n",len);
        return -1;
    }
    if (!(eph.sat=satno(SYS_GAL,prn))) {
        trace(2,"binex 0x01-14: satellite error prn=%d\n",prn);
        return -1;
    }
    if (!(eph.code&(1<<9))) return 0; /* only I/NAV */
    eph.A=sqrtA*sqrtA;
    eph.iodc=eph.iode;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,tocs);
    eph.ttr=adjweek(eph.toe,tow);
    if (ura < 0) {
        /*  -(SISA index + 1) for SISA index = 0-255 */
        eph.sva=-(ura + 1);
    }
    else {
        /* accuracy in distance in meters for SISA index = 0-125 */
        eph.sva=ura_index(ura,SYS_GAL);
    }
    
    if (!strstr(raw->opt,"-EPHALL")) {
        if (raw->nav.eph[eph.sat-1].iode==eph.iode&&
            raw->nav.eph[eph.sat-1].iodc==eph.iodc) return 0; /* unchanged */
    }
    raw->nav.eph[eph.sat-1]=eph;
    raw->ephsat=eph.sat;
    return 2;
}
/* decode binex mesaage 0x01: gnss navigaion informtion ----------------------*/
static int decode_bnx_01(raw_t *raw, unsigned char *buff, int len)
{
    char *msg;
    int srec=U1(buff),prn=U1(buff+1);
    
    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        prn=srec==0x01||srec==0x02||srec==0x04?prn+1:(srec==0x00?0:prn);
        sprintf(msg," subrec=%02X prn=%d",srec,prn);
    }
    switch (srec) {
        case 0x00: return decode_bnx_01_00(raw,buff+1,len-1);
        case 0x01: return decode_bnx_01_01(raw,buff+1,len-1);
        case 0x02: return decode_bnx_01_02(raw,buff+1,len-1);
        case 0x03: return decode_bnx_01_03(raw,buff+1,len-1);
        case 0x04: /* discard galileo ephemeris(binex) */
                   return -1;
        case 0x05: return decode_bnx_01_05(raw,buff+1,len-1);
        case 0x06: return decode_bnx_01_06(raw,buff+1,len-1);
        case 0x14: return decode_bnx_01_14(raw,buff+1,len-1);
    }
    return 0;
}
/* decode binex mesaage 0x02: generalized gnss data --------------------------*/
static int decode_bnx_02(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x02: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x03: generalized ancillary site data ----------------*/
static int decode_bnx_03(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x03: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7d: receiver internal state prototyping ------------*/
static int decode_bnx_7d(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x7d: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7e: ancillary site data prototyping ----------------*/
static int decode_bnx_7e(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x7e: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7f-00: jpl fiducial site ---------------------------*/
static int decode_bnx_7f_00(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x7f-00: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7f-01: ucar cosmic ---------------------------------*/
static int decode_bnx_7f_01(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x7f-01: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7f-02: trimble 4700 --------------------------------*/
static int decode_bnx_7f_02(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x7f-02: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7f-03: trimble netrs -------------------------------*/
static int decode_bnx_7f_03(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x7f-03: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7f-04: trimble netrs -------------------------------*/
static int decode_bnx_7f_04(raw_t *raw, unsigned char *buff, int len)
{
    trace(2,"binex 0x7f-04: not supported message\n");
    return 0;
}
/* decode binex mesaage 0x7f-05: trimble netr8 obs data ----------------------*/
static unsigned char *decode_bnx_7f_05_obs(raw_t *raw, unsigned char *buff,
                                           int sat, int nobs, obsd_t *data)
{
    const unsigned char codes_gps[32]={
        CODE_L1C ,CODE_L1C ,CODE_L1P ,CODE_L1W ,CODE_L1Y ,CODE_L1M , /*  0- 5 */
        CODE_L1X ,CODE_L1N ,CODE_NONE,CODE_NONE,CODE_L2W ,CODE_L2C , /*  6-11 */
        CODE_L2D ,CODE_L2S ,CODE_L2L ,CODE_L2X ,CODE_L2P ,CODE_L2W , /* 12-17 */
        CODE_L2Y ,CODE_L2M ,CODE_L2N ,CODE_NONE,CODE_NONE,CODE_L5X , /* 18-23 */
        CODE_L5I ,CODE_L5Q ,CODE_L5X                                 /* 24-26 */
    };
    const unsigned char codes_glo[32]={
        CODE_L1C ,CODE_L1C ,CODE_L1P ,CODE_NONE,CODE_NONE,CODE_NONE, /*  0- 5 */
        CODE_NONE,CODE_NONE,CODE_NONE,CODE_NONE,CODE_L2C ,CODE_L2C , /*  6-11 */
        CODE_L2P ,CODE_L3X ,CODE_L3I ,CODE_L3Q ,CODE_L3X             /* 12-16 */
    };
    const unsigned char codes_gal[32]={
        CODE_L1C ,CODE_L1A ,CODE_L1B ,CODE_L1C ,CODE_L1X ,CODE_L1Z , /*  0- 5 */
        CODE_L5X ,CODE_L5I ,CODE_L5Q ,CODE_L5X ,CODE_L7X ,CODE_L7I , /*  6-11 */
        CODE_L7Q ,CODE_L7X ,CODE_L8X ,CODE_L8I ,CODE_L8Q ,CODE_L8X , /* 12-17 */
        CODE_L6X ,CODE_L6A ,CODE_L6B ,CODE_L6C ,CODE_L6X ,CODE_L6Z , /* 18-23 */
    };
    const unsigned char codes_sbs[32]={
        CODE_L1C ,CODE_L1C ,CODE_NONE,CODE_NONE,CODE_NONE,CODE_NONE, /*  0- 5 */
        CODE_L5X ,CODE_L5I ,CODE_L5Q ,CODE_L5X                       /*  6- 9 */
    };
    const unsigned char codes_cmp[32]={
        CODE_L1X ,CODE_L1I ,CODE_L1Q ,CODE_L1X ,CODE_L7X ,CODE_L7I , /*  0- 5 */
        CODE_L7Q ,CODE_L7X ,CODE_L6X ,CODE_L6I ,CODE_L6Q ,CODE_L6X , /*  6-11 */
        CODE_L1X ,CODE_L1S ,CODE_L1L ,CODE_L1X                       /* 12-15 */
    };
    const unsigned char codes_qzs[32]={
        CODE_L1C ,CODE_L1C ,CODE_L1S ,CODE_L1L ,CODE_L1X ,CODE_NONE, /*  0- 5 */
        CODE_NONE,CODE_L2X ,CODE_L2S ,CODE_L2L ,CODE_L2X ,CODE_NONE, /*  6-11 */
        CODE_NONE,CODE_L5X ,CODE_L5I ,CODE_L5Q ,CODE_L5X ,CODE_NONE, /* 12-17 */
        CODE_NONE,CODE_L6X ,CODE_L6S ,CODE_L6L ,CODE_L6X ,CODE_NONE, /* 18-23 */
        CODE_NONE,CODE_NONE,CODE_NONE,CODE_NONE,CODE_NONE,CODE_NONE, /* 24-29 */
        CODE_L1Z                                                     /* 30-30 */
    };
    const unsigned char *codes=NULL;
    double range[8],phase[8],cnr[8],dopp[8]={0},acc,wl;
    unsigned char *p=buff;
    unsigned char flag,flags[4];
    int i,j,k,sys,fcn=-10,code[8],slip[8],pri[8],freq[8],slipcnt[8]={0},mask[8]={0};
    
    trace(5,"decode_bnx_7f_05_obs: sat=%2d nobs=%2d\n",sat,nobs);
    
    sys=satsys(sat,NULL);
    
    switch (sys) {
        case SYS_GPS: codes=codes_gps; break;
        case SYS_GLO: codes=codes_glo; break;
        case SYS_GAL: codes=codes_gal; break;
        case SYS_QZS: codes=codes_qzs; break;
        case SYS_SBS: codes=codes_sbs; break;
        case SYS_CMP: codes=codes_cmp; break;
    }
    for (i=0;i<nobs;i++) {
        
        flag   =getbitu(p,0,1);
        slip[i]=getbitu(p,2,1);
        code[i]=getbitu(p,3,5); p++;
        
        for (j=0;j<4;j++) flags[j]=0;
        
        for (j=0;flag&&j<4;j++) {
            flag=U1(p++);
            flags[flag&0x03]=flag&0x7F;
            flag&=0x80;
        }
        if (flags[2]) {
            fcn=getbits(flags+2,2,4);
        }
        acc=(flags[0]&0x20)?0.0001:0.00002; /* phase accuracy */
        
        cnr[i]=U1(p++)*0.4;
        
        if (i==0) {
            cnr[i]+=getbits(p,0,2)*0.1;
            range[i]=getbitu(p,2,32)*0.064+getbitu(p,34,6)*0.001; p+=5;
        }
        else if (flags[0]&0x40) {
            cnr[i]+=getbits(p,0,2)*0.1;
            range[i]=range[0]+getbits(p,4,20)*0.001; p+=3;
        }
        else {
            range[i]=range[0]+getbits(p,0,16)*0.001; p+=2;
        }
        if (flags[0]&0x40) {
            phase[i]=range[i]+getbits(p,0,24)*acc; p+=3;
        }
        else {
            cnr[i]+=getbits(p,0,2)*0.1;
            phase[i]=range[i]+getbits(p,2,22)*acc; p+=3;
        }
        if (flags[0]&0x04) {
            dopp[i]=getbits(p,0,24)/256.0; p+=3;
        }
        if (flags[0]&0x08) {
            if (flags[0]&0x10) {
                slipcnt[i]=U2(p); p+=2;
            }
            else {
                slipcnt[i]=U1(p); p+=1;
            }
        }
        trace(5,"(%d) CODE=%2d S=%d F=%02X %02X %02X %02X\n",i+1,
              code[i],slip[i],flags[0],flags[1],flags[2],flags[3]);
        trace(5,"(%d) P=%13.3f L=%13.3f D=%7.1f SNR=%4.1f SCNT=%2d\n",
              i+1,range[i],phase[i],dopp[i],cnr[i],slipcnt[i]);
    }
    if (!codes) {
        data->sat=0;
        return p;
    }
    data->time=raw->time;
    data->sat=sat;
    
    /* get code priority */
    for (i=0;i<nobs;i++) {
        code2obs(codes[code[i]&0x3F],freq+i);
        pri[i]=getcodepri(sys,codes[code[i]&0x3F],raw->opt);
        
        /* frequency index for beidou */
        if (sys==SYS_CMP) {
            if      (freq[i]==5) freq[i]=2; /* B2 */
            else if (freq[i]==4) freq[i]=3; /* B3 */
        }
    }
    for (i=0;i<NFREQ;i++) {
        for (j=0,k=-1;j<nobs;j++) {
            if (freq[j]==i+1&&(k<0||pri[j]>pri[k])) k=j;
        }
        if (k<0) {
            data->P[i]=data->L[i]=0.0;
            data->D[i]=0.0f;
            data->SNR[i]=data->LLI[i]=0;
            data->code[i]=CODE_NONE;
        }
        else {
            wl=satwavelen(sat,i,&raw->nav);
            if (sys==SYS_GLO&&fcn>=-7&&freq[k]<=2) {
                wl=CLIGHT/(freq[k]==1?FREQ1_GLO+DFRQ1_GLO*fcn:
                                      FREQ2_GLO+DFRQ2_GLO*fcn);
            }
            data->P[i]=range[k];
            data->L[i]=wl<=0.0?0.0:phase[k]/wl;
            data->D[i]=dopp[k];
            data->SNR[i]=(unsigned char)(cnr[k]/0.25+0.5);
            data->code[i]=codes[code[k]&0x3F];
            data->LLI[i]=slip[k]?1:0;
            mask[k]=1;
        }
    }
    for (;i<NFREQ+NEXOBS;i++) {
        for (k=0;k<nobs;k++) {
            if (!mask[k]) break;
        }
        if (k>=nobs) {
            data->P[i]=data->L[i]=0.0;
            data->D[i]=0.0f;
            data->SNR[i]=data->LLI[i]=0;
            data->code[i]=CODE_NONE;
        }
        else {
            wl=satwavelen(sat,freq[k]-1,&raw->nav);
            if (sys==SYS_GLO&&fcn>=-7&&freq[k]<=2) {
                wl=CLIGHT/(freq[k]==1?FREQ1_GLO+DFRQ1_GLO*fcn:
                                      FREQ2_GLO+DFRQ2_GLO*fcn);
            }
            data->P[i]=range[k];
            data->L[i]=wl<=0.0?0.0:phase[k]/wl;
            data->D[i]=dopp[k];
            data->SNR[i]=(unsigned char)(cnr[k]/0.25+0.5);
            data->code[i]=codes[code[k]&0x3F];
            data->LLI[i]=slip[k]?1:0;
            mask[k]=1;
        }
    }
    return p;
}
/* decode binex mesaage 0x7f-05: trimble netr8 -------------------------------*/
static int decode_bnx_7f_05(raw_t *raw, unsigned char *buff, int len)
{
    obsd_t data={{0}};
    double clkoff=0.0,toff[16]={0};
    char *msg;
    unsigned char *p=buff;
    unsigned int flag;
    int i,nsat,nobs,prn,sys,sat,clkrst=0,rsys=0,nsys=0,tsys[16]={0};
    
    trace(4,"decode_bnx_7f_05\n");
    
    raw->obs.n=0;
    flag=U1(p++);
    nsat=(int)(flag&0x3F)+1;
    
    if (flag&0x80) { /* rxclkoff */
        clkrst=getbitu(p,0, 2);
        clkoff=getbits(p,2,22)*1E-9; p+=3;
    }
    if (flag&0x40) { /* systime */
        nsys=getbitu(p,0,4);
        rsys=getbitu(p,4,4); p++;
        for (i=0;i<nsys;i++) {
            toff[i]=getbits(p,0,24)*1E-9;
            tsys[i]=getbitu(p,28,4); p+=4;
        }
    }
    for (i=0;i<nsat;i++) {
        prn =U1(p++);
        nobs=getbitu(p,1,3);
        sys =getbitu(p,4,4); p++;
        
        trace(5,"binex 0x7F-05 PRN=%3d SYS=%d NOBS=%d\n",prn,sys,nobs);
        
        switch (sys) {
            case 0: sat=satno(SYS_GPS,prn); break;
            case 1: sat=satno(SYS_GLO,prn); break;
            case 2: sat=satno(SYS_SBS,prn); break;
            case 3: sat=satno(SYS_GAL,prn); break;
            case 4: sat=satno(SYS_CMP,prn); break;
            case 5: sat=satno(SYS_QZS,prn); break;
            default: sat=0; break;
        }
        /* decode binex mesaage 0x7F-05 obs data */
        if (!(p=decode_bnx_7f_05_obs(raw,p,sat,nobs,&data))) return -1;
        
        if ((int)(p-buff)>len) {
            trace(2,"binex 0x7F-05 length error: nsat=%2d len=%d\n",nsat,len);
            return -1;
        }
        /* save obs data to obs buffer */
        if (data.sat&&raw->obs.n<MAXOBS) {
            raw->obs.data[raw->obs.n++]=data;
        }
    }
    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg," nsat=%2d",nsat);
    }
    return raw->obs.n>0?1:0;
}
/* decode binex mesaage 0x7f: gnss data prototyping --------------------------*/
static int decode_bnx_7f(raw_t *raw, unsigned char *buff, int len)
{
    const static double gpst0[]={1980,1,6,0,0,0};
    char *msg;
    unsigned char *p=buff;
    unsigned int srec,min,msec;
    
    srec=U1(p); p+=1; /* subrecord id */
    min =U4(p); p+=4;
    msec=U2(p); p+=2;
    raw->time=timeadd(epoch2time(gpst0),min*60.0+msec*0.001);
    
    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg," subrec=%02X time%s",srec,time_str(raw->time,3));
    }
    switch (srec) {
        case 0x00: return decode_bnx_7f_00(raw,buff+7,len-7);
        case 0x01: return decode_bnx_7f_01(raw,buff+7,len-7);
        case 0x02: return decode_bnx_7f_02(raw,buff+7,len-7);
        case 0x03: return decode_bnx_7f_03(raw,buff+7,len-7);
        case 0x04: return decode_bnx_7f_04(raw,buff+7,len-7);
        case 0x05: return decode_bnx_7f_05(raw,buff+7,len-7);
    }
    return 0;
}
/* decode binex mesaage ------------------------------------------------------*/
static int decode_bnx(raw_t *raw)
{
    unsigned int len,cs1,cs2;
    int rec,len_h;
    
    rec=raw->buff[1]; /* record id */
    
    /* record and header length */
    len_h=getbnxi(raw->buff+2,&len);
    
    trace(5,"decode_bnx: rec=%02x len=%d\n",rec,len);
    
    /* check parity */
    if (raw->len-1<128) {
        cs1=U1(raw->buff+raw->len);
        cs2=csum8(raw->buff+1,raw->len-1);
    }
    else {
        cs1=U2(raw->buff+raw->len);
        cs2=rtk_crc16(raw->buff+1,raw->len-1);
    }
    if (cs1!=cs2) {
        trace(2,"binex 0x%02X parity error CS=%X %X\n",rec,cs1,cs2);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype,"BINEX 0x%02X (%4d)",rec,raw->len);
    }
    /* decode binex message record */
    switch (rec) {
        case 0x00: return decode_bnx_00(raw,raw->buff+2+len_h,len);
        case 0x01: return decode_bnx_01(raw,raw->buff+2+len_h,len);
        case 0x02: return decode_bnx_02(raw,raw->buff+2+len_h,len);
        case 0x03: return decode_bnx_03(raw,raw->buff+2+len_h,len);
        case 0x7d: return decode_bnx_7d(raw,raw->buff+2+len_h,len);
        case 0x7e: return decode_bnx_7e(raw,raw->buff+2+len_h,len);
        case 0x7f: return decode_bnx_7f(raw,raw->buff+2+len_h,len);
    }
    return 0;
}
/* synchronize binex message -------------------------------------------------*/
static int sync_bnx(unsigned char *buff, unsigned char data)
{
    buff[0]=buff[1]; buff[1]=data;
    
    return buff[0]==BNXSYNC2&&
           (buff[1]==0x00||buff[1]==0x01||buff[1]==0x02||buff[1]==0x03||
            buff[1]==0x7D||buff[1]==0x7E||buff[1]==0x7F);
}
/* input binex message from stream ---------------------------------------------
* fetch next binex data and input a message from stream
* args   : raw_t *raw   IO     receiver raw data control struct
*          unsigned char data I stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris)
* notes  : support only the following message (ref [1])
*
*          - big-endian, regular CRC, forward record (sync=0xE2)
*          - record-subrecord:
*            0x01-01: decoded gps ephemeris
*            0x01-02: decoded glonass ephemeris
*            0x01-03: decoded sbas ephemeris
*            0x01-04: decoded galileo ephemeris
*            0x01-05: decoded beidou-2/compass ephemeris
*            0x01-06: decoded qzss ephemeris
*            0x7f-05: gnss data prototyping - trimble netr8
*
*          to specify input options, set rtcm->opt to the following option
*          strings separated by spaces.
*
*          -EPHALL  : input all ephemerides
*          -GLss    : select signal ss for GPS (ss=1C,1P,...)
*          -RLss    : select signal ss for GLO (ss=1C,1P,...)
*          -ELss    : select signal ss for GAL (ss=1C,1B,...)
*          -JLss    : select signal ss for QZS (ss=1C,2C,...)
*          -CLss    : select signal ss for BDS (ss=2I,2X,...)
*
*-----------------------------------------------------------------------------*/
extern int input_bnx(raw_t *raw, unsigned char data)
{
    unsigned int len;
    int len_h,len_c;
    
    trace(5,"input_bnx: data=%02x\n",data);
    
    /* synchronize binex message */
    if (raw->nbyte==0) {
        if (!sync_bnx(raw->buff,data)) return 0;
        raw->nbyte=2;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;
    if (raw->nbyte<4) return 0;
    
    len_h=getbnxi(raw->buff+2,&len);
    
    raw->len=len+len_h+2; /* length without crc */
    
    if (raw->len-1>4096) {
        trace(2,"binex length error: len=%d\n",raw->len-1);
        raw->nbyte=0;
        return -1;
    }
    len_c=raw->len-1<128?1:2;
    
    if (raw->nbyte<(int)(raw->len+len_c)) return 0;
    raw->nbyte=0;
    
    /* decode binex message */
    return decode_bnx(raw);
}
/* input binex message from file -----------------------------------------------
* fetch next binex data and input a message from file
* args   : raw_t  *raw   IO     receiver raw data control struct
*          FILE   *fp    I      file pointer
* return : status(-2: end of file, -1...9: same as above)
*-----------------------------------------------------------------------------*/
extern int input_bnxf(raw_t *raw, FILE *fp)
{
    unsigned int len;
    int i,data,len_h,len_c;
    
    trace(4,"input_bnxf\n");
    
    if (raw->nbyte==0) {
        for (i=0;;i++) {
            if ((data=fgetc(fp))==EOF) return -2;
            if (sync_bnx(raw->buff,(unsigned char)data)) break;
            if (i>=4096) return 0;
        }
    }
    if (fread(raw->buff+2,1,4,fp)<4) return -2;
    
    len_h=getbnxi(raw->buff+2,&len);
    
    raw->len=len+len_h+2;
    
    if (raw->len-1>4096) {
        trace(2,"binex length error: len=%d\n",raw->len-1);
        raw->nbyte=0;
        return -1;
    }
    len_c=raw->len-1<128?1:2;
    
    if (fread(raw->buff+6,1,raw->len+len_c-6,fp)<(size_t)(raw->len+len_c-6)) {
        return -2;
    }
    raw->nbyte=0;
    
    /* decode binex message */
    return decode_bnx(raw);
}

#define SQR(x)      ((x)*(x))
#define NUMSYS      6                   /* number of systems */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */
#define MINFREQ_GLO -7                  /* min frequency number glonass */
#define MAXFREQ_GLO 13                  /* max frequency number glonass */
#define NINCOBS     262144              /* inclimental number of obs data */
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    int frq[MAXOBSTYPE];                /* signal frequency (1:L1,2:L2,...) */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    unsigned char pri [MAXOBSTYPE];     /* signal priority (15-0) */
    unsigned char type[MAXOBSTYPE];     /* type (0:C,1:L,2:D,3:S) */
    unsigned char code[MAXOBSTYPE];     /* obs code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;
static void setstr(char *dst, const char *src, int n)
{
    char *p=dst;
    const char *q=src;
    while (*q&&q<src+n) *p++=*q++;
    *p--='\0';
    while (p>=dst&&*p==' ') *p--='\0';
}
static void init_sta(sta_t *sta)
{
    int i;
    *sta->name   ='\0';
    *sta->marker ='\0';
    *sta->antdes ='\0';
    *sta->antsno ='\0';
    *sta->rectype='\0';
    *sta->recver ='\0';
    *sta->recsno ='\0';
    sta->antsetup=sta->itrf=sta->deltype=0;
    for (i=0;i<3;i++) sta->pos[i]=0.0;
    for (i=0;i<3;i++) sta->del[i]=0.0;
    sta->hgt=0.0;
}
static void saveslips(unsigned char slips[][NFREQ], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ;i++) {
        if (data->LLI[i]&1) slips[data->sat-1][i]|=1;
    }
}
/* restore slips -------------------------------------------------------------*/
static void restslips(unsigned char slips[][NFREQ], obsd_t *data)
{
    int i;
    for (i=0;i<NFREQ;i++) {
        if (slips[data->sat-1][i]&1) data->LLI[i]|=1;
        slips[data->sat-1][i]=0;
    }
}
/* add obs data --------------------------------------------------------------*/
static int addobsdata(obs_t *obs, const obsd_t *data)
{
    obsd_t *obs_data;
    
    if (obs->nmax<=obs->n) {
        if (obs->nmax<=0) obs->nmax=NINCOBS; else obs->nmax*=2;
        if (!(obs_data=(obsd_t *)realloc(obs->data,sizeof(obsd_t)*obs->nmax))) {
            trace(1,"addobsdata: memalloc error n=%dx%d\n",sizeof(obsd_t),obs->nmax);
            free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
            return -1;
        }
        obs->data=obs_data;
    }
    obs->data[obs->n++]=*data;
    return 1;
}
typedef struct {                /* stream file type */
    int    format;              /* stream format (STRFMT_???) */
    int    sat;                 /* input satellite */
    obs_t  *obs;                /* input observation data */
    nav_t  *nav;                /* input navigation data */
    gtime_t time;               /* current time */
    rtcm_t rtcm;                /* rtcm data */
    raw_t  raw;                 /* receiver raw data */
    rnxctr_t rnx;               /* rinex data */
    FILE   *fp;                 /* file pointer */
} strfile_t;
static int add_eph(nav_t *nav, const eph_t *eph)
{
    eph_t *nav_eph;
    
    if (nav->nmax<=nav->n) {
        nav->nmax+=1024;
        if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->nmax))) {
            trace(1,"decode_eph malloc error: n=%d\n",nav->nmax);
            free(nav->eph); nav->eph=NULL; nav->n=nav->nmax=0;
            return 0;
        }
        nav->eph=nav_eph;
    }
    nav->eph[nav->n++]=*eph;
    return 1;
}
static int add_geph(nav_t *nav, const geph_t *geph)
{
    geph_t *nav_geph;
    
    if (nav->ngmax<=nav->ng) {
        nav->ngmax+=1024;
        if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ngmax))) {
            trace(1,"decode_geph malloc error: n=%d\n",nav->ngmax);
            free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
            return 0;
        }
        nav->geph=nav_geph;
    }
    nav->geph[nav->ng++]=*geph;
    return 1;
}
static int add_seph(nav_t *nav, const seph_t *seph)
{
    seph_t *nav_seph;
    
    if (nav->nsmax<=nav->ns) {
        nav->nsmax+=1024;
        if (!(nav_seph=(seph_t *)realloc(nav->seph,sizeof(seph_t)*nav->nsmax))) {
            trace(1,"decode_seph malloc error: n=%d\n",nav->nsmax);
            free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
            return 0;
        }
        nav->seph=nav_seph;
    }
    nav->seph[nav->ns++]=*seph;
    return 1;
}
/* set system mask -----------------------------------------------------------*/
static int set_sysmask(const char *opt)
{
    const char *p;
    int mask=SYS_NONE;
    
    if (!(p=strstr(opt,"-SYS="))) return SYS_ALL;
    
    for (p+=5;*p&&*p!=' ';p++) {
        switch (*p) {
            case 'G': mask|=SYS_GPS; break;
            case 'R': mask|=SYS_GLO; break;
            case 'E': mask|=SYS_GAL; break;
            case 'J': mask|=SYS_QZS; break;
            case 'C': mask|=SYS_CMP; break;
            case 'S': mask|=SYS_SBS; break;
        }
    }
    return mask;
}
static int readbnxfp(FILE *fp, gtime_t ts, gtime_t te, double tint,
                     const char *opt, int flag, int index, char *type,
                     obs_t *obs, nav_t *nav, sta_t *sta, const prcopt_t *popt)
{
    unsigned char slips[MAXSAT][NFREQ]={{0}};
    int stat;
    static raw_t raw; /* static declaration to reduction stack */
    unsigned int len;
    int data,len_h,len_c;
    int dec_type;
    int sys, sat, prn;
    int mask;
    int i;
    
    memset(&raw, '\0', sizeof(raw));
    
    /* initialize receiver raw and rtcm control */
    init_raw(&raw);
    strcpy(raw.opt,opt);
    /* set system mask */
    /* mask=set_sysmask(opt); */
    mask = popt->navsys;
    
    while (1) {
        if ((data=fgetc(fp))==EOF) break;
        if (!sync_bnx(raw.buff,(unsigned char)data)) continue;
        
        if (fread(raw.buff+2,1,4,fp)<4) break;
        
        len_h=getbnxi(raw.buff+2,&len);
        
        raw.len=len+len_h+2;
        
        if (raw.len-1>4096) {
            trace(2,"binex length error: len=%d\n",raw.len-1);
            raw.nbyte=0;
            continue;
        }
        len_c=raw.len-1<128?1:2;
        
        if (raw.len+len_c-6<1) {
            trace(2,"binex length error: len=%d\n",raw.len-1);
            raw.nbyte=0;
            continue;
        }
        
        if (fread(raw.buff+6,1,raw.len+len_c-6,fp)<(size_t)(raw.len+len_c-6)) {
            break;
        }
        raw.nbyte=0;
        if ((dec_type=decode_bnx(&raw))>=1) {
            /* avioid duplicated if overlapped data */
#if 0
            if (tend.time&&timediff(str->time,tend)<=0.0) continue;
#endif

            switch (dec_type) {
                case  1:
                    for (i=0;i<raw.obs.n;i++) {
                        /* save cycle-slip */
                        saveslips(slips,raw.obs.data+i);
                    }
                    /* screen data by time */
                    if (raw.obs.n>0&&!screent(raw.obs.data[0].time,ts,te,tint)) continue;
                    
                    for (i=0;i<raw.obs.n;i++) {
                        /* restore cycle-slip */
                        restslips(slips,raw.obs.data+i);
                        
                        raw.obs.data[i].rcv=(unsigned char)index;
                        
                        if (!(satsys(raw.obs.data[i].sat,NULL)&mask)) {
                            continue;
                        }
                        
                        /* save obs data */
                        if ((stat=addobsdata(obs,raw.obs.data+i))<0) break;
                    }
                    break;
                case  2:
                    sat = raw.ephsat;
                    
                    sys = satsys(sat,&prn);
                    
                    if (!(sys&mask)) {
                        break;
                    }
                    
                    switch (sys) {
                        case SYS_GLO:
                            stat = add_geph(nav,&raw.nav.geph[prn-1]);
                            break;
                        case SYS_SBS:
                            stat = add_seph(nav,&raw.nav.seph[prn-MINPRNSBS]);
                            break;
                        default:
                            stat = add_eph(nav,&raw.nav.eph[sat-1]);
                            break;
                    }
                    if (!stat) return 0;
                    break;
#if 0
                case  3: convsbs(ofp,opt,str,n);       break;
                case 31: convlex(ofp,opt,str,n);       break;
#endif
            }
#if 0
            te=str->time; if (ts.time==0) ts=te;
            
            /* set approx position */
            if (type==1&&!opt->autopos&&norm(opt->apppos,3)<=0.0) {
                setapppos(str,opt);
            }
            if (opt->te.time&&timediff(te,opt->te)>10.0) break;
#endif
        }
    }
    
    free_raw(&raw);
    
    return 1;
}
static int readbnxfile(const char *file, gtime_t ts, gtime_t te, double tint,
                       const char *opt, int flag, int index, char *type,
                       obs_t *obs, nav_t *nav, sta_t *sta, const prcopt_t *popt)
{
    FILE *fp;
    int cstat,stat;
    char tmpfile[1024];
    
    trace(3,"readbnxfile: file=%s flag=%d index=%d\n",file,flag,index);
    
    if (sta) init_sta(sta);
    
    /* uncompress file */
    if ((cstat=uncompress(file,tmpfile))<0) {
        trace(2,"binex file uncompact error: %s\n",file);
        return 0;
    }
    if (!(fp=fopen(cstat?tmpfile:file,"rb"))) {
        trace(2,"binex file open error: %s\n",cstat?tmpfile:file);
        return 0;
    }
    /* read binex file */
    stat=readbnxfp(fp,ts,te,tint,opt,flag,index,type,obs,nav,sta,popt);
    
    fclose(fp);
    
    /* delete temporary file */
    if (cstat) remove(tmpfile);
    
    return stat;
}
extern int readbnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
                    double tint, const char *opt, obs_t *obs, nav_t *nav,
                    sta_t *sta, const prcopt_t *popt)
{
    int i,n,stat=0;
    const char *p;
    char type=' ',*files[MAXEXFILE]={0};
    
    trace(3,"readbnxt: file=%s rcv=%d\n",file,rcv);
    
    if (!*file) {
        return readbnxfp(stdin,ts,te,tint,opt,0,1,&type,obs,nav,sta,popt);
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]);
            return -1;
        }
    }
    /* expand wild-card */
    if ((n=expath(file,files,MAXEXFILE))<=0) {
        for (i=0;i<MAXEXFILE;i++) free(files[i]);
        return 0;
    }
    /* read binex files */
    for (i=0;i<n&&stat>=0;i++) {
        stat=readbnxfile(files[i],ts,te,tint,opt,0,rcv,&type,obs,nav,sta,popt);
    }
    /* if station name empty, set 4-char name from file head */
    if (sta) {
        if (!(p=strrchr(file,FILEPATHSEP))) p=file-1;
        if (!*sta->name) setstr(sta->name,p+1,4);
    }
    for (i=0;i<MAXEXFILE;i++) free(files[i]);
    
    return stat;
}
