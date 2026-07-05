# Makefile — FalconSign SamplerZ simulation with Verilator.
#
#   make test   — generate memory, compile, simulate
#   make mem    — (re)generate memory.hex only
#   make sim    — compile + run against existing memory.hex
#   make clean

VERILATOR ?= verilator
PYTHON    ?= python3
INC       := ../include
SAM       := ../sampler

SRCS := \
    $(INC)/time_counter.sv \
    $(INC)/fpr_mul_opt_2pipe.sv \
    $(INC)/fpr_mul_const_opt_2pipe.sv \
    $(INC)/floating_point_add_2pip.sv \
    $(INC)/floating_point_sub_2pip.sv \
    $(INC)/fp_sub_s.sv \
    $(INC)/fp_mult_s.sv \
    $(INC)/fp_flt2i_int32_s.sv \
    $(INC)/fp_i2flt_int32_s.sv \
    $(INC)/fp_i2flt_int32.sv \
    $(INC)/fp_flt2i_int64_s.sv \
    $(INC)/fp_i2flt_int64_s.sv \
    $(SAM)/chacha20.sv \
    $(SAM)/berexp.sv \
    $(SAM)/fpr_cal.sv \
    $(SAM)/pre_samp.sv \
    $(SAM)/samp_loop.sv \
    $(SAM)/refill_control.sv \
    $(SAM)/samplerz.sv \
    $(SAM)/samplerz_top.sv \
    samplerz_tb.sv

VFLAGS := --binary --timing --sv \
          -I$(INC) \
          -Wno-WIDTHTRUNC -Wno-WIDTHEXPAND -Wno-TIMESCALEMOD \
          -Wno-MULTITOP -Wno-INITIALDLY -Wno-CASEINCOMPLETE \
          -Wno-UNOPTFLAT -Wno-UNSIGNED -Wno-LITENDIAN \
          -Wno-WIDTHXZEXPAND \
          --top-module samplerz_tb -o Vsamplerz_tb

STAT_SRCS := $(filter-out samplerz_tb.sv, $(SRCS)) samplerz_stat_tb.sv

STAT_VFLAGS := $(VFLAGS) --top-module samplerz_stat_tb -o Vsamplerz_stat_tb

.PHONY: test mem sim stat_sim stat clean

test: mem sim

mem:
	$(PYTHON) gen_mem.py memory.hex

sim: $(SRCS)
	$(VERILATOR) $(VFLAGS) $(SRCS)
	./obj_dir/Vsamplerz_tb

stat_sim: $(STAT_SRCS)
	$(VERILATOR) $(STAT_VFLAGS) $(STAT_SRCS)

stat: mem stat_sim
	$(PYTHON) stat_test.py

clean:
	rm -rf obj_dir memory.hex
