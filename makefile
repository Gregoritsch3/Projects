CC = g++

CFLAGS = -c -Wall

all: prog

prog: RadioactiveDecay.o newtonraphson.o
	$(CC) -o prog RadioactiveDecay.o newtonraphson.o

RadioactiveDecay.o: RadioactiveDecay.cpp
	$(CC) $(CFLAGS) RadioactiveDecay.cpp

newtonraphson.o: newtonraphson.cpp
	$(CC) $(CFLAGS) newtonraphson.cpp

clean:
	rm -rf *.o
