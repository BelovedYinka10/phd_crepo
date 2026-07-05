# Makefile — build and run the FP module tests with Verilator.
#
#   make test     — generate vectors then simulate
#   make vectors  — (re)generate fp_vectors.txt only
#   make sim      — compile + run testbench against existing vectors
#   make clean

VERILATOR ?= verilator
PYTHON    ?= python3
NVEC      ?= 500
INC       := ../include

# All RTL sources needed (order matters: dependencies before dependents)
SRCS := \
    $(INC)/fpr_mul_opt_2pipe.sv \
    $(INC)/floating_point_add_2pip.sv \
    $(INC)/floating_point_sub_2pip.sv \
    $(INC)/fp_sub_s.sv \
    $(INC)/fp_mult_s.sv \
    $(INC)/fp_flt2i_int32_s.sv \
    $(INC)/fp_i2flt_int32_s.sv \
    $(INC)/fp_i2flt_int32.sv \
    $(INC)/fp_flt2i_int64_s.sv \
    $(INC)/fp_i2flt_int64_s.sv \
    fp_tb.sv

VFLAGS := --binary --timing --sv \
          -Wno-WIDTHTRUNC -Wno-WIDTHEXPAND -Wno-INITIALDLY \
          -Wno-CASEINCOMPLETE -Wno-UNOPTFLAT \
          --top-module fp_tb -o Vfp_tb

.PHONY: test vectors sim clean

test: vectors sim

vectors:
	$(PYTHON) gen_fp_vectors.py $(NVEC)

sim: $(SRCS)
	$(VERILATOR) $(VFLAGS) $(SRCS)
	./obj_dir/Vfp_tb

clean:
	rm -rf obj_dir fp_vectors.txt
