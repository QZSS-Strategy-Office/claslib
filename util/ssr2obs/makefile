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

CONFOBS	 = sample.conf
TIME1	 = -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59 -ti 1
NAV1	 = ../../data/20190827.nav
CSSR1	 = ../../data/20190827q.l6
VRS1	 = vrs2019239Q.obs -r
TIME2	 = -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59 -ti 1
NAV2	 = ../../data/20191215.nav
CSSR2	 = ../../data/20191215b.l6
VRS2	 = vrs20191215b.obs -r
RTCM2	 = vrs20191215b.rtcm3 -b
test: ssr2obs testobs1 testobs2 testrtcm2

testdump: # Dump of L6 message
	./ssr2obs -k $(CONFOBS) -dump $(CSSR2)

testobs1: # Output in RINEX 3 format
	./ssr2obs -k $(CONFOBS) $(TIME1) $(NAV1) $(CSSR1) -o $(VRS1)

testobs2: # Output in RINEX 3 format(ST12)
	./ssr2obs -k $(CONFOBS) $(TIME2) $(NAV2) $(CSSR2) -o $(VRS2)

testrtcm2: # Output in RTCM 3 MSM format
	./ssr2obs -k $(CONFOBS) $(TIME2) $(NAV2) $(CSSR2) -o $(RTCM2)
