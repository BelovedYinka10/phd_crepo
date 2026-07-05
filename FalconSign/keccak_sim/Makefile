# Makefile — SHAKE256 / Keccak-f[1600] Vivado XSim simulation
#
# Prerequisites:
#   source /tools/Xilinx/Vivado/<version>/settings64.sh
#
# Targets:
#   make sim      — generate vectors, compile, elaborate, run (batch)
#   make project  — create .xpr Vivado GUI project
#   make vectors  — (re)generate shake256_vectors.hex only
#   make clean

VIVADO  ?= vivado
PYTHON  ?= python3

ROOT    := ..
KECCAK  := $(ROOT)/keccak
SIM     := .

SRCS := \
    $(KECCAK)/keccak_f1600.sv \
    $(KECCAK)/shake256.sv \
    $(SIM)/shake256_tb.sv

TOP     := shake256_tb
SNAP    := shake256_tb_snap

XVLOG_FLAGS := -sv --nolog -Wno-WFA_INCOMPATIBLE_TYPES
XELAB_FLAGS := -debug typical --nolog -Wno-CASEINCOMPLETE
XSIM_FLAGS  := --runall --nolog

.PHONY: sim project vectors clean

sim: vectors compile elaborate run

vectors:
	$(PYTHON) gen_vectors.py

compile: $(SRCS)
	@echo "=== xvlog: compiling ==="
	xvlog $(XVLOG_FLAGS) $(SRCS)

elaborate:
	@echo "=== xelab: top=$(TOP) ==="
	xelab $(XELAB_FLAGS) $(TOP) -s $(SNAP)

run:
	@echo "=== xsim: running ==="
	xsim $(SNAP) $(XSIM_FLAGS)

project: vectors
	$(VIVADO) -mode batch -source create_project.tcl \
	    -log create_project.log -nojournal

clean:
	rm -rf xsim.dir .Xil *.log *.jou *.pb *.wdb *.str
	rm -f shake256_vectors.hex
	rm -rf falcon_shake256
