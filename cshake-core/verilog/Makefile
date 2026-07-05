# Makefile — cshake-core Verilator simulation (macOS / Linux)
#
# Step 1: generate test vectors (run from testvectors_generator/):
#   make vectors TARGET=cshake256   (or shake256, cshake128, shake128)
#
# Step 2: simulate:
#   make sim
#   make clean

PYTHON  ?= python3
CXX     ?= g++
NUM     ?= 10

SRCS := clog2.v keccak_pkg.v rc.v stateram_inference.v state_ram.v \
        transform.v keccak_top.v control_path.v data_path.v \
        tb.v testvectors.v

VL_FLAGS := --binary --timing \
    -Wno-WIDTHTRUNC -Wno-WIDTHEXPAND -Wno-CASEINCOMPLETE \
    -Wno-UNOPTFLAT -Wno-ASCRANGE -Wno-LATCH

.PHONY: sim vectors clean

sim: obj_dir/Vcshake_tb
	./obj_dir/Vcshake_tb

obj_dir/Vcshake_tb: $(SRCS)
	verilator $(VL_FLAGS) --top-module cshake_tb -o Vcshake_tb $(SRCS)

# Regenerate testvectors.v (TARGET = shake256 | shake128 | cshake256 | cshake128)
TARGET ?= cshake256
vectors:
	cd ../testvectors_generator && \
	    $(CXX) -D$(shell echo $(TARGET) | tr a-z A-Z) -DNUMTESTS=$(NUM) -DVERILOG \
	        generator.cpp fips202.c -o generator && \
	    ./generator > ../verilog/testvectors.v
	@echo "testvectors.v written for $(TARGET) with NUMTESTS=$(NUM)"

clean:
	rm -rf obj_dir ../testvectors_generator/generator
