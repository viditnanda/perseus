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

# the discrete-gradient cubical engine against the classic one (barcodes over F_2, F_3, F_5, Q; certificates)
test-cubical: tests/t_cubical_morse.cpp $(HDR)
	$(CXX) -std=c++17 -O2 -Wall -pthread $(INC) tests/t_cubical_morse.cpp -o t_cubical_morse
	./t_cubical_morse 300
	PERSEUS_LONG_COLUMN=0 ./t_cubical_morse 100
	$(CXX) -std=c++17 -O2 -Wall $(INC) tests/t_bitcol.cpp -o t_bitcol
	./t_bitcol 100

clean:
	rm -f perseus t_cubical_morse t_bitcol
