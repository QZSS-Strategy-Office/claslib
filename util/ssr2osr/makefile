# makefile for ssr2osr

CC		 = gcc
CFLAGS	 = \
			-DENAGAL \
			-DENAQZS \
			-DENA_SSR2OSR \
			-DLAPACK \
			-DNFREQ=3 \
			-DTRACE \
			-I$(DIRCLAS) \
			-I$(DIRRTKLIB) \
			-O3 \
			-Wall \
			-g \

LIBS	 = \
			-lblas \
			-llapack \
			-lm \

RTKLIB	 = \
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

GSILIB	 = \
			isb.c \

CLAS	 = \
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
TARGETS	 = ssr2osr

all: $(TARGETS)

ssr2osr: $(OBJS) $(DIROBJ)/ssr2osr.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

ssr2obs: $(OBJS) $(DIROBJ)/ssr2obs.o
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
	rm -rf $(DIROBJ) $(TARGETS) *.csv *.llh*

.PHONY: all clean test testdump testnmea1 testnmea2

CONFOSR	 = sample.conf
TIME1	 = -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59 -ti 1
TIME1	 =
OBS1	 = ../../data/20190827q.obs
NAV1	 = ../../data/20190827.nav
CSSR1	 = ../../data/20190827q.l6
LLH1	 = $(basename $(notdir $(OBS1))).llh
TIME2	 = -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59 -ti 1
TIME2	 =
OBS2	 = ../../data/20191215b.bnx
CSSR2	 = ../../data/20191215b.l6
LLH2	 = $(basename $(notdir $(OBS2))).llh

test: ssr2osr testnmea1 testnmea2

testdump : # Dump of L6 message
	./ssr2osr -k $(CONFOSR) -dump $(CSSR2)

testnmea1 : # Output corrections represented in observation space (OSR)
	./ssr2osr -k $(CONFOSR) $(TIME1) $(OBS1) $(NAV1) $(CSSR1) -o $(LLH1)

testnmea2 : # Output corrections represented in observation space (OSR-ST12)
	./ssr2osr -k $(CONFOSR) $(TIME2) $(OBS2) $(CSSR2) -o $(LLH2)
