#OPT=-O3 -Wall -DSIMPSON -DLAPLACIAN2
OPT=-O3 -Wall

all: sample1.x sample1b.x sample1c.x sample1d.x sample2a.x sample2b.x sample2c.x sample3.x sample4.x sample5.x sample6.x sample7.x sample8.x

sample1.x: sample1.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample1b.x: sample1b.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample1c.x: sample1c.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample1d.x: sample1d.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample2a.x: sample2a.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample2b.x: sample2b.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample2c.x: sample2c.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample3.x: sample3.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample4.x: sample4.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample5.x: sample5.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample6.x: sample6.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample7.x: sample7.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample8.x: sample8.cc simplebounce.o
	g++ $^ -o $@ $(OPT)

simplebounce.o: simplebounce.cc
	g++ $^ -c $(OPT)

clean:
	rm *.x *.o -f

