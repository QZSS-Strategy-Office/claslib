# makefile for ssr2obs

CC		 = gcc
CFLAGS	 = \
			-DCSSR2OSR_VRS \
			-DENAGAL \
			-DENAQZS \
			-DLAPACK \
			-DNEXOBS=2 \
			-DNFREQ=3 \
			-DTRACE \
			-I$(DIRCLAS) \
			-I$(DIRRTKLIB) \
			-O3 \
			-Wall \
			-g \

LIBS   = \
			-lblas \
			-llapack \
			-lm \
			-lpthread \

RTKLIB	  = \
			ephemeris.c \
			geoid.c \
			ionex.c \
			lambda.c \
			options.c \
			pntpos.c \
			postpos.c \
			ppp.c \
			ppp_ar.c \
			preceph.c \
			qzslex.c \
			rcvraw.c \
			rinex.c \
			rtcm.c \
			rtcm2.c \
			rtcm3.c \
			rtcm3e.c \
			rtkcmn.c \
			rtkpos.c \
			sbas.c \
			solution.c \

GSILIB	  = \
			isb.c \

CLAS	  = \
			cssr.c \
			cssr2osr.c \
			grid.c \
			stec.c \

RCV		 = \
			binex.c \
			crescent.c \
			gw10.c \
			javad.c \
			novatel.c \
			nvs.c \
			rcvlex.c \
			rt17.c \
			skytraq.c \
			ss2.c \
			ublox.c \

DIRCLAS	 = ../../src/claslib
DIRGSILIB= ../../src/gsilib
DIRRCV	 = ../../src/rcv
DIRRTKLIB= ../../src/rtklib
DIROBJ	 = obj
OBJS	 = $(addprefix $(DIROBJ)/, \
			$(CLAS:.c=.o) \
			$(GSILIB:.c=.o) \
			$(RCV:.c=.o) \
			$(RTKLIB:.c=.o))
TARGETS	 = ssr2obs

all: $(TARGETS)

ssr2obs: $(OBJS) $(DIROBJ)/ssr2obs.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

ssr2osr: $(OBJS) $(DIROBJ)/ssr2osr.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(DIROBJ)/%.o: $(DIRCLAS)/%.c
	mkdir -p $(DIROBJ)
	$(CC) -c -o $@ $< $(CFLAGS)

$(DIROBJ)/%.o: $(DIRRCV)/%.c
	mkdir -p $(DIROBJ)
	$(CC) -c -o $@ $< $(CFLAGS)

$(DIROBJ)/%.o: $(DIRRTKLIB)/%.c
	mkdir -p $(DIROBJ)
	$(CC) -c -o $@ $< $(CFLAGS)

$(DIROBJ)/%.o: $(DIRGSILIB)/%.c
	mkdir -p $(DIROBJ)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -rf $(DIROBJ) $(TARGETS) *.csv *.obs *.rtcm3

.PHONY: all clean test testdump testobs1 testobs2 testrtcm2

CONF    = sample.conf
CSSR1	= ../../data/2019239Q.l6
CSSR2	= ../../data/2019349B.l6
NAV1	= ../../data/sept_2019239.nav
NAV2	= ../../data/sept_2019349.nav
OUTPUT1	= $(basename $(notdir $(CSSR1))).obs   -r
OUTPUT2	= $(basename $(notdir $(CSSR2))).obs   -r
OUTPUT3	= $(basename $(notdir $(CSSR2))).rtcm3 -b
TIME1	= -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59 -ti 1
TIME2	= -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59 -ti 1
test: ssr2obs testobs1 testobs2 testrtcm2
testdump: # Dump of L6 message
	./ssr2obs -k $(CONF) -dump $(CSSR2)
testobs1: # Output in RINEX 3 format
	./ssr2obs -k $(CONF) $(NAV1) $(CSSR1) -o $(OUTPUT1) $(TIME1)
testobs2: # Output in RINEX 3 format(ST12)
	./ssr2obs -k $(CONF) $(NAV2) $(CSSR2) -o $(OUTPUT2) $(TIME2)
testrtcm2: # Output in RTCM 3 MSM format
	./ssr2obs -k $(CONF) $(NAV2) $(CSSR2) -o $(OUTPUT3) $(TIME2)

