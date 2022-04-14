
all: main 
build: NWC_math.o BitOperate.o NWC.o

CC = g++
INSTDIR = /usr/local/bin
INCLUDE = .
CFLAGS = -g -std=c++11 -Wall -ansi
LIBS += -framework CoreFoundation

main: main.o main_CF.o build		
	$(CC) -o main.exe main.o NWC_math.o	BitOperate.o NWC.o
	$(CC) -o main_CF.exe main_CF.o NWC_math.o	BitOperate.o NWC.o

main.o: main.cpp
main_CF.o: main_CF.cpp
NWC_math.o: NWC_math.cpp NWC_math.h
BitOperate.o: BitOperate.cpp BitOperate.h
NWC.o: NWC.cpp NWC.h

clean:
	rm *.exe
	rm *.o
	rm *.txt

