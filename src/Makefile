G++ = g++
G++_FLAGS = -c -std=c++11 -Wall

OBJECTS = GaussianPyramid.o OpticalFlow.o Stochastic.o
TARGET = flow.a

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar rcs $(TARGET) $(OBJECTS) $(LD_FLAGS)

%.o : %.cpp
	$(G++) $(G++_FLAGS) $<

clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all clean

