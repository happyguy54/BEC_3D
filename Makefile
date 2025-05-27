CC=g++
CFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags`
LIBS=`root-config --libs`

DEPS=fit.h c2.h utils.h

all: fit clean

%.o: %.cxx $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

fit: fit.o c2.o utils.o
	$(CC) -g -o fit fit.o c2.o utils.o $(LDFLAGS) $(LIBS) -lMinuit

clean:
	rm *.o
