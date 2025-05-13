CXX = g++
CXXFLAGS = -std=c++20 -O2

TARGET = simpleks

all: $(TARGET)

simpleks : main.o simpleks.o
	$(CXX) -o $@ $^

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

simpleks.o: simpleks.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f main main.o simpleks.o

.PHONY: all clean
