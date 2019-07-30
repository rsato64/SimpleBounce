OPT=-O3 -Wall 

all: sample1.x sample2c.x sample7.x sample8.x

sample1.x: sample1.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample2c.x: sample2c.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample7.x: sample7.cc simplebounce.o
	g++ $^ -o $@ $(OPT)
sample8.x: sample8.cc simplebounce.o
	g++ $^ -o $@ $(OPT)

simplebounce.o: simplebounce.cc
	g++ $^ -c $(OPT)

clean:
	rm *.x *.o -f

