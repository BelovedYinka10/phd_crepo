# Makefile — FalconSign SamplerZ: Vivado XSim simulation (Vivado 2019+)
#
# Prerequisites:
#   - Source Vivado settings first:
#       source /tools/Xilinx/Vivado/2023.2/settings64.sh
#   - Python 3 for memory generation
#
# Targets:
#   make sim      — compile + elaborate + run simulation (batch, no GUI)
#   make project  — create .xpr project (open in Vivado GUI, click Run Sim)
#   make mem      — (re)generate memory.hex
#   make wave     — open existing waveform in Vivado GUI
#   make clean

VIVADO  ?= vivado
PYTHON  ?= python3

# Relative to this Makefile (FalconSign/vivado/)
ROOT    := ..
INC     := $(ROOT)/include
SAM     := $(ROOT)/sampler
SIM     := $(ROOT)/sampler_sim

# All source files — ORDER MATTERS for xvlog (packages first, tb last)
SRCS := \
    $(INC)/falconsoar_pkg.sv \
    $(INC)/sample_pkg.sv \
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
    $(INC)/counter_ce.sv \
    $(INC)/rbsh.sv \
    $(SAM)/chacha20.sv \
    $(SAM)/berexp.sv \
    $(SAM)/fpr_cal.sv \
    $(SAM)/pre_samp.sv \
    $(SAM)/samp_loop.sv \
    $(SAM)/refill_control.sv \
    $(SAM)/samplerz.sv \
    $(SAM)/samplerz_top.sv \
    $(SIM)/samplerz_tb.sv

TOP     := samplerz_tb
SNAP    := samplerz_tb_snap

# xvlog flags:
#   -sv           = SystemVerilog mode
#   -i $(INC)     = include search path for `include directives
#   --nolog       = suppress xvlog.log (we write our own)
#   -d XSIM_SIM   = optional define so RTL can distinguish simulation
XVLOG_FLAGS := -sv -i $(INC) --nolog -d XSIM_SIM \
               -Wno-WFA_INCOMPATIBLE_TYPES

XELAB_FLAGS := -debug typical --nolog \
               -Wno-CASEINCOMPLETE

XSIM_FLAGS  := --runall --nolog

.PHONY: sim project wave mem clean

# ── Batch simulation ───────────────────────────────────────────────────────────
sim: mem compile elaborate run

mem:
	$(PYTHON) $(SIM)/gen_mem.py memory.hex

compile: $(SRCS)
	@echo "=== xvlog: compiling $(words $(SRCS)) source files ==="
	xvlog $(XVLOG_FLAGS) $(SRCS)

elaborate:
	@echo "=== xelab: elaborating top=$(TOP) ==="
	xelab $(XELAB_FLAGS) $(TOP) -s $(SNAP)

run:
	@echo "=== xsim: running simulation ==="
	xsim $(SNAP) $(XSIM_FLAGS)

# ── Create Vivado GUI project ─────────────────────────────────────────────────
project: mem
	$(VIVADO) -mode batch -source create_project.tcl \
	    -log create_project.log -nojournal

# ── Open existing waveform in Vivado GUI ──────────────────────────────────────
wave:
	$(VIVADO) -mode gui -source open_wave.tcl &

clean:
	rm -rf xsim.dir .Xil *.log *.jou *.pb *.wdb *.str memory.hex
	rm -rf falcon_sampler
