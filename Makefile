all : read_ems
CXXFLAGS=`root-config --cflags` -I./ -I$(HOME)/work/ -std=c++14
LIBS= `root-config --libs` -lMinuit -lboost_program_options -lfmt 

read_ems: read_ems.o
	g++ $(LIBS) read_ems.o $(LIBS) -o read_ems


read_ems.o :  read_ems.cpp multistep.h
	g++ -c read_ems.cpp $(CXXFLAGS)  -o read_ems.o

#.o.cpp:
#	g++  -o $@ $(CXXFLAGS) -c $^
	
clean :
				rm -f *.o rm *.so* *.a read_ems
