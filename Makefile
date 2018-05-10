GXX:= g++
FLAGS:= -O -std=c++11

all: test

test: main.o primes.o l.a
	$(GXX) main.o primes.o l.a -o test

primes.o: primes.cxx
	$(GXX) $(FLAGS) -c primes.cxx

main.o: main.cpp
	$(GXX) $(FLAGS) -c main.cpp

clear:
	rm -f *.o test
