CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Iinclude
LDFLAGS ?= -lm

SRC := src/image.cpp src/fusion.cpp src/main.cpp
BIN := build/exposure_fusion

all: $(BIN)

$(BIN): $(SRC) include/exposure_fusion.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(SRC) -o $(BIN) $(LDFLAGS)

run: $(BIN)
	./$(BIN) --inputs data/inputs/belgium/BelgiumHouse_A.png,data/inputs/belgium/BelgiumHouse_B.png --outdir data/output/belgium --levels 8

clean:
	rm -rf build

.PHONY: all run clean
