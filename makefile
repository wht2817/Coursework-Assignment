CXX = mpicxx
CXXFLAGS = -Wall -O3
HDRS = CommandLineOptions.h SPH.h
OBJS = main.o CommandLineOptions.o SPH.o
LDLIBS = -lblas -lboost_program_options
TARGET = myprog

default: $(TARGET)
all: $(TARGET)

%.o : %.cpp $(HDRS)
	$(CXX) -g $(CXXFLAGS) -o $@ -c $<
	
$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDLIBS)

.PHONY: clean

clean:
	rm -f $(TARGET) *.o


