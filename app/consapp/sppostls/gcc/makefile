# makefile for sppostls
SRC     = ../../../../src

OPTS    = -DTRACE -DENAGLO -DENAQZS -DENAGAL -DENACMP -DENAIRN -DNFREQ=5 -DNEXOBS=2
CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTS) -g
LDLIBS  = ../../../../lib/iers/gcc/iers.a -lm

all        : sppostls

sppostls   : sppostls.o postls.o rtkcmn.o lambda.o rinex.o pntpos.o ephemeris.o
sppostls   : preceph.o sbas.o ionex.o geoid.o solution.o options.o tides.o gpt3.o

sppostls.o : ../sppostls.c
	$(CC) -c $(CFLAGS) ../sppostls.c
postls.o   : $(SRC)/postls.c
	$(CC) -c $(CFLAGS) $(SRC)/postls.c
rtkcmn.o   : $(SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(SRC)/rtkcmn.c
lambda.o   : $(SRC)/lambda.c
	$(CC) -c $(CFLAGS) $(SRC)/lambda.c
rinex.o    : $(SRC)/rinex.c
	$(CC) -c $(CFLAGS) $(SRC)/rinex.c
pntpos.o   : $(SRC)/pntpos.c
	$(CC) -c $(CFLAGS) $(SRC)/pntpos.c
ephemeris.o: $(SRC)/ephemeris.c
	$(CC) -c $(CFLAGS) $(SRC)/ephemeris.c
preceph.o  : $(SRC)/preceph.c
	$(CC) -c $(CFLAGS) $(SRC)/preceph.c
sbas.o     : $(SRC)/sbas.c
	$(CC) -c $(CFLAGS) $(SRC)/sbas.c
ionex.o    : $(SRC)/ionex.c
	$(CC) -c $(CFLAGS) $(SRC)/ionex.c
geoid.o    : $(SRC)/geoid.c
	$(CC) -c $(CFLAGS) $(SRC)/geoid.c
solution.o : $(SRC)/solution.c
	$(CC) -c $(CFLAGS) $(SRC)/solution.c
options.o  : $(SRC)/options.c
	$(CC) -c $(CFLAGS) $(SRC)/options.c
tides.o    : $(SRC)/tides.c
	$(CC) -c $(CFLAGS) $(SRC)/tides.c
gpt3.o     : $(SRC)/gpt3.c
	$(CC) -c $(CFLAGS) $(SRC)/gpt3.c

sppostls.o postls.o rtkcmn.o lambda.o rinex.o pntpos.o ephemeris.o preceph.o sbas.o ionex.o geoid.o solution.o options.o tides.o gpt3.o : $(SRC)/rtklib.h

clean :
	rm -f sppostls sppostls.exe *.o *.pos *.trace
