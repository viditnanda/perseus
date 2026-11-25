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

# the discrete-gradient and planar cubical engines against the classic one (barcodes over F_2, F_3, F_5, Q;
# certificates)
test-cubical: tests/t_cubical_morse.cpp tests/t_cubical_planar.cpp $(HDR)
	$(CXX) -std=c++17 -O2 -Wall -pthread $(INC) tests/t_cubical_morse.cpp -o t_cubical_morse
	./t_cubical_morse 300
	PERSEUS_LONG_COLUMN=0 ./t_cubical_morse 100
	$(CXX) -std=c++17 -O2 -Wall $(INC) tests/t_bitcol.cpp -o t_bitcol
	./t_bitcol 100
	$(CXX) -std=c++17 -O2 -Wall -pthread $(INC) tests/t_cubical_planar.cpp -o t_cubical_planar
	./t_cubical_planar 400
	PERSEUS_PLANAR_COMPRESS=1 ./t_cubical_planar 200 2
	$(CXX) -std=c++17 -O2 -Wall $(INC) tests/t_cubes.cpp -o t_cubes
	./t_cubes 600

# chain complexes over Z, simplicial sets, integral outputs, magnitude homology, Laurent coefficients
test-chain: tests/t_chain.cpp tests/t_magnitude.cpp tests/t_laurent.cpp tests/t_howell.cpp $(HDR)
	$(CXX) -std=c++17 -O2 -Wall -pthread $(INC) tests/t_chain.cpp -o t_chain
	$(CXX) -std=c++17 -O2 -Wall -pthread $(INC) tests/t_magnitude.cpp -o t_magnitude
	$(CXX) -std=c++17 -O2 -Wall -pthread $(INC) tests/t_laurent.cpp -o t_laurent
	./t_chain 200 && ./t_magnitude && ./t_laurent 60
	$(CXX) -std=c++17 -O2 -Wall $(INC) tests/t_exact.cpp -o t_exact
	./t_exact 2000
	$(CXX) -std=c++17 -O2 -Wall -Wextra $(INC) tests/t_howell.cpp -o t_howell
	./t_howell 2000

test-rips-components: tests/t_rips_components.cpp $(HDR)
	$(CXX) -std=c++17 -O2 -pthread $(INC) $< -o t_rips_components
	./t_rips_components

# finite solvable groups: Anick resolution, integral homology, mod-p cohomology rings (closed forms, Kunneth,
# universal coefficients, HAP and Green-King ground truth)
test-pcgroup: tests/t_pcgroup.cpp $(HDR)
	$(CXX) -std=c++17 -O2 -Wall $(INC) tests/t_pcgroup.cpp -o t_pcgroup
	./t_pcgroup

clean:
	rm -f perseus t_rips_components t_cubical_morse t_bitcol t_cubical_planar t_chain t_magnitude t_laurent t_exact t_howell t_pcgroup
