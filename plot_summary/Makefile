CC=g++
CFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags`
LIBS=`root-config --libs`

DEPS=plot.h utils.h

all: plot clean

%.o: %.cxx $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

plot: plot.o utils.o
	$(CC) -g -o plot plot.o utils.o $(LDFLAGS) $(LIBS)

clean:
	rm *.o
