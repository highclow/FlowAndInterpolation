GOOGLE_TEST_LIB = gtest
GOOGLE_TEST_INCLUDE = /usr/include

G++ = g++
G++_FLAGS = -c -std=c++11 -Wall -I $(GOOGLE_TEST_INCLUDE) -I ../src/
LD_FLAGS = -L /usr/lib -l $(GOOGLE_TEST_LIB) -l pthread -L ../src/ -l flow

OBJECTS = gtest_main.o test_image.o
TARGET = gtest

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(G++) -o $(TARGET) $(OBJECTS) $(LD_FLAGS)

%.o : %.cpp
	$(G++) $(G++_FLAGS) $<

clean:
	rm -f $(TARGET) $(OBJECTS)
                    
.PHONY: all clean
