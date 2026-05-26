# makefile for rnx2rtkp

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
OPTS   += -DENA_PPP_RTK
OPTS   += -DENA_REL_VRS




SRCS  = rtkcmn.c rinex.c rtkpos.c postpos.c solution.c lambda.c geoid.c sbas.c \
		preceph.c pntpos.c ephemeris.c options.c ppp.c ppp_ar.c rtcm.c rtcm2.c \
		rtcm3.c rtcm3e.c ionex.c qzslex.c rcvraw.c
SRCS += stec.c
SRCS += isb.c
SRCS += grid.c cssr.c
SRCS += ppprtk.c cssr2osr.c
SRCS += rtkvrs.c

RCV_SRCS = binex.c gw10.c novatel.c skytraq.c ublox.c crescent.c javad.c nvs.c \
		   ss2.c rt17.c

OBJS = $(SRCS:.c=.o)
RCV_OBJS = $(RCV_SRCS:.c=.o)

all        : rnx2rtkp

-include makefile.depend

rnx2rtkp.o : rnx2rtkp.c
	$(CC) -c $(CFLAGS) $<

rnx2rtkp   : rnx2rtkp.o $(OBJS) $(RCV_OBJS)

%.o : $(SRC)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

%.o : $(RCV_SRC)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

CMD = ./rnx2rtkp$(EXT)

# for kinematic
CONF_L6   = kinematic.conf
CONF_VRS  = kinematic_vrs.conf


# for static
#CONF_L6  = static.conf
#CONF_VRS = static_vrs.conf


OPTT   = -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59
OBS    = $(DATADIR)0627239Q.obs
VRS    = $(DATADIR)vrs2019239Q.obs
NAV    = $(DATADIR)sept_2019239.nav
L6     = $(DATADIR)2019239Q.l6
OPTT_W = -ts 2018/11/25 0:00:00 -te 2018/11/25 00:59:59 -l6w 2028
OBS_W  = $(DATADIR)0161329A.obs
NAV_W  = $(DATADIR)tskc2018329.nav
L6_W   = $(DATADIR)2018328X_329A.l6
OPTT_B = -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59
OBS_B  = $(DATADIR)0627239Q.bnx
NAV_B  = $(DATADIR)sept_2019239.nav
L6_B   = $(DATADIR)2019239Q.l6
OPTT_ST12 = -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59
OBS_ST12_B= $(DATADIR)0627349AB.bnx
NAV_ST12  = $(DATADIR)sept_2019349.nav
L6_ST12   = $(DATADIR)2019349B.l6
VRS_ST12= $(DATADIR)vrs2019349B.obs

INPUT_L6   = $(OBS) $(NAV) $(L6)
INPUT_VRS  = $(OBS) $(NAV) $(VRS)
OUTPUT_L6  = $(basename $(notdir $(OBS)))_L6.nmea
OUTPUT_VRS = $(basename $(notdir $(OBS)))_VRS.nmea
INPUT_L6W  = $(OBS_W) $(NAV_W) $(L6_W)
INPUT_L6B  = $(OBS_B) $(NAV_B) $(L6_B)
INPUT_ST12 = $(OBS_ST12_B) $(NAV_ST12) $(L6_ST12)
INPUT_NONAV= $(OBS_ST12_B) $(L6_ST12)
INPUT_VRS_ST12 = $(OBS_ST12_B) $(NAV_ST12) $(VRS_ST12)
OUTPUT_L6  = $(basename $(notdir $(OBS)))_claslib_L6.nmea
OUTPUT_VRS = $(basename $(notdir $(OBS)))_claslib_VRS.nmea
OUTPUT_L6W = $(basename $(notdir $(OBS_W)))_claslib_L6W.nmea
OUTPUT_L6B = $(basename $(notdir $(OBS_B)))_claslib_L6B.nmea
OUTPUT_ST12= $(basename $(notdir $(OBS_ST12_B)))_claslib_ST12.nmea
OUTPUT_NONAV= $(basename $(notdir $(OBS_ST12_B)))_claslib_NONAV.nmea
OUTPUT_VRS_ST12 = $(basename $(notdir $(OBS_ST12_B)))_claslib_VRS_ST12.nmea

test_L6 :
	$(CMD) -ti 1 $(OPTT) -x 2 -k $(CONF_L6) $(INPUT_L6) -o $(OUTPUT_L6)

test_L6_week:
	$(CMD) -ti 1 $(OPTT_W) -x 2 -k $(CONF_L6) $(INPUT_L6W) -o $(OUTPUT_L6W)

test_L6_bnx:
	$(CMD) -ti 1 $(OPTT_B) -x 2 -k $(CONF_L6) $(INPUT_L6B) -o $(OUTPUT_L6B)

test_VRS :
	$(CMD) -ti 1 $(OPTT) -x 2 -k $(CONF_VRS) $(INPUT_VRS) -o $(OUTPUT_VRS)

test_ST12:
	$(CMD) -ti 1 $(OPTT_ST12) -x 2 -k $(CONF_L6) $(INPUT_ST12) -o $(OUTPUT_ST12)
	
test_NONAV:
	$(CMD) -ti 1 $(OPTT_ST12) -x 2 -k $(CONF_L6) $(INPUT_NONAV) -o $(OUTPUT_NONAV)

test_ST12_VRS:
	$(CMD) -ti 1 $(OPTT_ST12) -x 3 -k $(CONF_VRS) $(INPUT_VRS_ST12) -o $(OUTPUT_VRS_ST12)

clean :
	rm -f rnx2rtkp$(EXT) *.o

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
