CC = g++
CFLAGS = -std=c++0x -I/opt/boost/boost_1_57_0 -Wall -Wextra -pedantic -g3 -DDEBUG
LFLAGS = -L/opt/boost/boost_1_57_0/stage/lib -lboost_serialization -O3

#all:map_ser

OBJS=Driver.o Solution.o

rod: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LFLAGS)

%.o : %.cpp
	$(CC) -c $(CFLAGS) $^

#map_ser.o:map_ser.cpp
#	$(CC) -c $(CFLAGS) map_ser.cpp


clean:
	rm *.o rod
