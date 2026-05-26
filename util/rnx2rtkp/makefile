# makefile for rnx2rtkp

CC		 = gcc
CFLAGS	 = \
			-DENA_PPP_RTK \
			-DENA_REL_VRS \
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
			ppprtk.c \
			rtkvrs.c \
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
TARGETS	 = rnx2rtkp

all: $(TARGETS)

ssr2obs: $(OBJS) $(DIROBJ)/ssr2obs.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

ssr2osr: $(OBJS) $(DIROBJ)/ssr2osr.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

rnx2rtkp: $(OBJS) $(DIROBJ)/rnx2rtkp.o
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

.PHONY: all clean test test_L6 test_L6_week test_L6_bnx test_VRS test_ST12 test_NONAV test_ST12_VRS

CONF1       = kinematic.conf
CONF2       = kinematic_vrs.conf
CSSR1       = ../../data/2018328X_329A.l6
CSSR2       = ../../data/2019239Q.l6
CSSR3       = ../../data/2019349B.l6
NAV1        = ../../data/tskc2018329.nav
NAV2        = ../../data/sept_2019239.nav
NAV3        = ../../data/sept_2019349.nav
OBS1        = ../../data/0161329A.obs
OBS2        = ../../data/0627239Q.obs
OBS2B       = ../../data/0627239Q.bnx
OBS3        = ../../data/0627349AB.bnx
TIME1       = -ti 1 -ts 2018/11/25 00:00:00 -te 2018/11/25 00:59:59 -l6w 2028
TIME2       = -ti 1 -ts 2019/08/27 16:00:00 -te 2019/08/27 16:59:59
TIME3       = -ti 1 -ts 2019/12/15 01:00:00 -te 2019/12/15 01:59:59
VRS1        = ../../data/vrs2019239Q.obs
VRS2        = ../../data/vrs2019349B.obs
INPUT_L6W   = $(OBS1)  $(NAV1) $(CSSR1)
INPUT_L6    = $(OBS2)  $(NAV2) $(CSSR2)
INPUT_L6B   = $(OBS2B) $(NAV2) $(CSSR2)
INPUT_NONAV = $(OBS3)  $(CSSR3)
INPUT_ST12  = $(OBS3)  $(NAV3) $(CSSR3)
INPUT_VRS1  = $(OBS2)  $(NAV2) $(VRS1)
INPUT_VRS2  = $(OBS3)  $(NAV3) $(VRS2)
OUTPUT_L6W  = $(basename $(notdir $(OBS1)))_L6W.pos
OUTPUT_L6   = $(basename $(notdir $(OBS2)))_L6.pos
OUTPUT_L6B  = $(basename $(notdir $(OBS2B)))_L6B.pos
OUTPUT_NONAV= $(basename $(notdir $(OBS3)))_NONAV.pos
OUTPUT_ST12 = $(basename $(notdir $(OBS3)))_ST12.pos
OUTPUT_VRS1 = $(basename $(notdir $(OBS2)))_VRS1.pos
OUTPUT_VRS2 = $(basename $(notdir $(OBS3)))_VRS2.pos
test: test_L6 test_L6_week test_L6_bnx test_VRS test_ST12 test_NONAV test_ST12_VRS
test_L6 :
	./rnx2rtkp -k $(CONF1) $(INPUT_L6)    -o $(OUTPUT_L6)    # $(TIME2)
test_L6_week:
	./rnx2rtkp -k $(CONF1) $(INPUT_L6W)   -o $(OUTPUT_L6W)   # $(TIME1)
test_L6_bnx:
	./rnx2rtkp -k $(CONF1) $(INPUT_L6B)   -o $(OUTPUT_L6B)   # $(TIME2)
test_VRS :
	./rnx2rtkp -k $(CONF2) $(INPUT_VRS1)  -o $(OUTPUT_VRS1)  # $(TIME2)
test_ST12:
	./rnx2rtkp -k $(CONF1) $(INPUT_ST12)  -o $(OUTPUT_ST12)  # $(TIME3)
test_NONAV:
	./rnx2rtkp -k $(CONF1) $(INPUT_NONAV) -o $(OUTPUT_NONAV) # $(TIME3)
test_ST12_VRS:
	./rnx2rtkp -k $(CONF2) $(INPUT_VRS2)  -o $(OUTPUT_VRS2)  # $(TIME3)

