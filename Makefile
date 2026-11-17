CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra -Wno-unused-parameter -pthread
INC = -Iinclude
HDR = $(wildcard include/perseus/*.hpp)

all: perseus

perseus: src/perseus.cpp $(HDR)
	$(CXX) $(CXXFLAGS) $(INC) src/perseus.cpp -o perseus

portable: src/perseus.cpp $(HDR)
	$(CXX) $(subst -march=native,,$(CXXFLAGS)) $(INC) src/perseus.cpp -o perseus

test: perseus
	python3 tests/run_tests.py

clean:
	rm -f perseus
