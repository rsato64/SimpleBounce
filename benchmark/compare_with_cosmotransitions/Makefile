include ../../Makefile.inc

all: model1_190103714.x model2_190103714.x model3_190103714.x model4_190103714.x model5_190103714.x model6_190103714.x model7_190103714.x model8_190103714.x model1a_190610829.x model1b_190610829.x model2a_190610829.x model2b_190610829.x

%.x: %.cc ../../simplebounce.o
	$(CXX) $^ -o $@ $(OPT) -I../../

clean:
	rm *.x *.o -f

