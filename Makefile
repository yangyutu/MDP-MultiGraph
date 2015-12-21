CC = g++
CFLAGS = -std=c++0x -I/opt/boost/boost_1_57_0 -Wall -Wextra -pedantic -O3 -DDEBUG
LFLAGS = -L/opt/boost/boost_1_57_0/stage/lib -lboost_serialization

#all:map_ser

OBJS=Driver.o Solution.o

rod: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LFLAGS)

%.o : %.cpp
	$(CC) -c $(CFLAGS) $^

#map_ser.o:map_ser.cpp
#	$(CC) -c $(CFLAGS) map_ser.cpp

rod_static: $(OBJS)
	$(CC) -o $@ $(OBJS) -static $(LFLAGS)

clean:
	rm *.o rod
