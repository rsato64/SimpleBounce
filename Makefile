include Makefile.inc

all: sample1.x
	make -C benchmark

%.x: %.cc simplebounce.o
	$(CXX) $^ -o $@ $(OPT)

simplebounce.o: simplebounce.cc
	$(CXX) $^ -c $(OPT)

clean:
	rm *.x *.o -f
	make -C benchmark clean

