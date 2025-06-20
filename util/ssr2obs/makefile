# makefile for ssr2obs

ifeq ($(OS),Windows_NT)
    # for MinGW
    CC = i686-w64-mingw32-gcc
    EXT=.exe
    CFLAGS  = -Wall -O3 -ansi -pedantic -I$(SRC) $(OPTS)
    CFLAGS += -DWIN32
    LDLIBS  = -lm
    LDLIBS += -lwinmm
    DATADIR = ..\\..\\data\\
else
    # for gcc
    CC = gcc
    EXT=
    # without lapack
    #CFLAGS  = -Wall -O3 -ansi -pedantic -I$(SRC) $(OPTS)
    #LDLIBS  = -lm
    # with lapack
    CFLAGS  = -Wall -O3 -ansi -pedantic -I$(SRC) $(OPTS)
    CFLAGS += -DLAPACK
    LDLIBS  = -lm
    LDLIBS += -llapack -lblas
    DATADIR = ../../data/
endif


SRC     = ../../src
RCV_SRC = ../../src/rcv

OPTS    = -DTRACE -DENAGAL -DENAQZS -DNFREQ=3
OPTS   += -DNEXOBS=2 -DDLL
OPTS   += -DCSSR2OSR_VRS



SRCS  = rtkcmn.c rinex.c rtkpos.c postpos.c solution.c lambda.c geoid.c sbas.c \
		preceph.c pntpos.c ephemeris.c options.c ppp.c ppp_ar.c rtcm.c rtcm2.c \
		rtcm3.c rtcm3e.c ionex.c qzslex.c rcvraw.c
SRCS += stec.c
SRCS += isb.c
SRCS += grid.c cssr.c
SRCS += cssr2osr.c

RCV_SRCS = binex.c novatel.c ublox.c

OBJS = $(SRCS:.c=.o)
RCV_OBJS = $(RCV_SRCS:.c=.o)

EXE = ssr2obs

all : $(EXE)

-include makefile.depend

ssr2obs.o : ssr2obs.c
	$(CC) -c $(CFLAGS) $<

ssr2obs   : ssr2obs.o $(OBJS) $(RCV_OBJS)

%.o : $(SRC)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

%.o : $(RCV_SRC)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

CMD    = ./ssr2obs$(EXT)

NAV    = $(DATADIR)sept_2019239.nav
CSSR   = $(DATADIR)2019239Q.l6
NAV_ST12 = $(DATADIR)sept_2019349.nav
CSSR_ST12 = $(DATADIR)2019349B.l6
L6_CH1  = $(DATADIR)2025157U_ch1.l6
L6_CH2  = $(DATADIR)2025157U_ch2.l6
NAV_2CH = $(DATADIR)0627157U.nav

# Dump of L6 message 
test1d:
	$(CMD) -k sample.conf -dump $(CSSR_ST12)

# Output in RINEX 3 format
test1r:
	$(CMD) -k sample.conf -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59 -ti 1 $(NAV) $(CSSR) -r -o $(DATADIR)vrs2019239Q.obs

# Output in RINEX 3 format(ST12)
test1r_ST12:
	$(CMD) -k sample.conf -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59 -ti 1 $(NAV_ST12) $(CSSR_ST12) -r -o $(DATADIR)vrs2019349B.obs

# Output in RINEX 3 format(CH1)
test1r_CH1:
	$(CMD) -k sample.conf -ts 2025/06/6 20:00:00 -te 2025/06/6 20:59:59 -ti 1 $(NAV_2CH) $(L6_CH1) -r -o $(DATADIR)2025157U_ch1.obs

# Output in RINEX 3 format(CH2)
test1r_CH2:
	$(CMD) -k sample.conf -ts 2025/06/6 20:00:00 -te 2025/06/6 20:59:59 -ti 1 $(NAV_2CH) $(L6_CH2) -r -o $(DATADIR)2025157U_ch2.obs

# Output in RTCM 3 MSM format
test1m:
	$(CMD) -k sample.conf -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59 -ti 1 $(NAV_ST12) $(CSSR_ST12) -b -o $(DATADIR)vrs2019349B.rtcm3


clean :
	rm -f ssr2obs$(EXT) *.o

cleano :
	rm -f *.o

depend :
	-@ rm -f makefile.depend
	-@ for i in $(SRCS:%.c=$(SRC)/%.c); do\
		$(CC) -MM -MG $(CFLAGS) $$i >> makefile.depend;\
	done
	-@ for i in $(RCV_SRCS:%.c=$(RCV_SRC)/%.c); do\
		$(CC) -MM -MG $(CFLAGS) $$i >> makefile.depend;\
	done
