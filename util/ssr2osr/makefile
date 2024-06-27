# makefile for ssr2osr

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
OPTS   += -DENA_SSR2OSR


SRCS  = rtkcmn.c rinex.c rtkpos.c postpos.c solution.c lambda.c geoid.c sbas.c \
		preceph.c pntpos.c ephemeris.c options.c ppp.c ppp_ar.c rtcm.c rtcm2.c \
		rtcm3.c rtcm3e.c ionex.c qzslex.c rcvraw.c
SRCS += stec.c
SRCS += isb.c
SRCS += grid.c cssr.c
SRCS += cssr2osr.c ppprtk.c

RCV_SRCS = binex.c gw10.c novatel.c skytraq.c ublox.c crescent.c javad.c nvs.c \
		   ss2.c rt17.c

OBJS = $(SRCS:.c=.o)
RCV_OBJS = $(RCV_SRCS:.c=.o)

all        : ssr2osr

-include makefile.depend

ssr2osr.o : ssr2osr.c
	$(CC) -c $(CFLAGS) $<

ssr2osr   : ssr2osr.o $(OBJS) $(RCV_OBJS)

%.o : $(SRC)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

%.o : $(RCV_SRC)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

CMD = ./ssr2osr$(EXT)

CONF   = sample.conf

OPTT   = -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59
OBS    = $(DATADIR)0627239Q.obs
NAV    = $(DATADIR)sept_2019239.nav
CSSR   = $(DATADIR)2019239Q.l6
OPTT_ST12 = -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59
OBS_ST12  = $(DATADIR)0627349AB.bnx
NAV_ST12  = $(DATADIR)sept_2019349.nav
CSSR_ST12 = $(DATADIR)2019349B.l6

INPUT  = $(OBS) $(NAV) $(CSSR)
OUTPUT = $(basename $(notdir $(OBS))).nmea
INPUT_ST12  = $(OBS_ST12) $(NAV_ST12) $(CSSR_ST12)
OUTPUT_ST12 = $(basename $(notdir $(OBS_ST12)))_ST12.nmea

# Dump of L6 message 
test1d :
	$(CMD) -k sample.conf -dump $(CSSR_ST12)
	
# Output corrections represented in observation space (OSR)
test1 :
	$(CMD) -ti 1 $(OPTT) -x 2 -k $(CONF) $(INPUT) -o $(OUTPUT)
# Output corrections represented in observation space (OSR) (ST12)
test1_ST12 :
	$(CMD) -ti 1 $(OPTT_ST12) -x 2 -k $(CONF) $(INPUT_ST12) -o $(OUTPUT_ST12)



clean :
	rm -f ssr2osr$(EXT) *.o

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
