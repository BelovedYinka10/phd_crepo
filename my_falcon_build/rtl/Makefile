# Makefile — build & verify the Falcon BaseSampler RTL with Verilator.
#
#   make test     # generate vectors from the Python reference, then simulate
#   make vectors  # (re)generate vectors.txt only
#   make sim      # build + run the testbench against existing vectors.txt
#   make clean

VERILATOR ?= verilator
PYTHON    ?= python3
TOP        = base_sampler_tb
SRCS       = base_sampler.v base_sampler_tb.v
NVEC      ?= 2000

.PHONY: test vectors sim clean

test: vectors sim

vectors:
	$(PYTHON) gen_vectors.py $(NVEC)

sim: $(SRCS)
	$(VERILATOR) --binary --timing -Wno-WIDTHTRUNC -Wno-WIDTHEXPAND \
	    --top-module $(TOP) -o V$(TOP) $(SRCS)
	./obj_dir/V$(TOP)

clean:
	rm -rf obj_dir vectors.txt
