include ../../Makefile.inc

all: change_n_tau.x

%.x: %.cc ../../simplebounce.o
	$(CXX) $^ -o $@ $(OPT) -I../../

clean:
	rm *.x *.o *.dat *.pdf -f

