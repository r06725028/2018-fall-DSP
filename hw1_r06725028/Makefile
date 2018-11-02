.PHONY: all clean

CFLAGS+=
LDFLAGS+=-lm     # link to math library

#TARGET=test

#all: $(TARGET)
# type make/make all to compile test_hmm
all:
	gcc train.cpp -o train
	gcc test.cpp -o test

run:
	./train 1000 model_init.txt seq_model_01.txt model_01.txt
	./train 1000 model_init.txt seq_model_02.txt model_02.txt
	./train 1000 model_init.txt seq_model_03.txt model_03.txt
	./train 1000 model_init.txt seq_model_04.txt model_04.txt
	./train 1000 model_init.txt seq_model_05.txt model_05.txt
	./test modellist.txt testing_data1.txt result1.txt
	./test modellist.txt testing_data2.txt result2.txt

clean:
	rm ./train
	rm ./test
	#$(RM) $(TARGET)   # type make clean to remove the compiled file
