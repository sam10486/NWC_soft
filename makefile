
all: main 
build: NWC_math.o BitOperate.o NWC.o

CC = g++
INSTDIR = /usr/local/bin
INCLUDE = .
CFLAGS = -g -std=c++11 -Wall -ansi
LIBS += -framework CoreFoundation

main: main.o main_CF.o main_AE.o main_BU.o main_R16_BU.o main_barrett.o main_radix_test.o \
		test_file.o test_file_mem_AE.o check_AE_NWC.o build		
	$(CC) -o main.exe main.o NWC_math.o	BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o main_CF.exe main_CF.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o main_AE.exe main_AE.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o main_BU.exe main_BU.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o main_R16_BU.exe main_R16_BU.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o main_barrett.exe main_barrett.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o main_radix_test.exe main_radix_test.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o test_file.exe test_file.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o test_file_mem_AE.exe test_file_mem_AE.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm
	$(CC) -o check_AE_NWC.exe check_AE_NWC.o NWC_math.o BitOperate.o NWC.o -lntl -lgmp -lm


main.o: main.cpp
main_CF.o: main_CF.cpp
main_AE.o: main_AE.cpp
main_BU.o: main_BU.cpp
main_R16_BU.o: main_R16_BU.cpp
main_barrett.o: main_barrett.cpp
main_radix_test.o: main_radix_test.cpp
test_file.o: test_file.cpp
test_file_mem_AE.o: test_file_mem_AE.cpp
check_AE_NWC.o: check_AE_NWC.cpp

NWC_math.o: NWC_math.cpp NWC_math.h
BitOperate.o: BitOperate.cpp BitOperate.h
NWC.o: NWC.cpp NWC.h

clean:
	rm *.exe
	rm *.o
	rm *.txt

