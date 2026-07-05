# create_project.tcl — Vivado GUI project for SHAKE256 / Keccak-f[1600]
#
# Run once:
#   cd FalconSign/keccak_sim
#   python3 gen_vectors.py
#   vivado -mode batch -source create_project.tcl
#
# Then open falcon_shake256/falcon_shake256.xpr and click:
#   Flow Navigator → Simulation → Run Simulation → Run Behavioral Simulation

set SCRIPT_DIR [file normalize [file dirname [info script]]]
set ROOT       [file normalize "$SCRIPT_DIR/.."]
set KECCAK     "$ROOT/keccak"
set SIM        "$SCRIPT_DIR"
set PROJ_DIR   "$SCRIPT_DIR/falcon_shake256"

# ── Generate vectors if missing ───────────────────────────────────────────────
if {![file exists "$SIM/shake256_vectors.hex"]} {
    puts "Generating shake256_vectors.hex..."
    exec python3 "$SIM/gen_vectors.py"
}

# ── Create project ────────────────────────────────────────────────────────────
# Part: ZCU104 (Zynq UltraScale+). Change PART to match your board.
create_project falcon_shake256 $PROJ_DIR \
    -part xczu7ev-ffvc1156-2-e -force

set_property default_lib       work          [current_project]
set_property target_language   SystemVerilog [current_project]
set_property simulator_language Mixed        [current_project]

# ── Design sources ────────────────────────────────────────────────────────────
add_files -norecurse [list \
    $KECCAK/keccak_f1600.sv \
    $KECCAK/shake256.sv     \
]

# ── Testbench ─────────────────────────────────────────────────────────────────
add_files -fileset sim_1 -norecurse $SIM/shake256_tb.sv

# Copy vectors file where XSim can find it
file copy -force "$SIM/shake256_vectors.hex" "$SCRIPT_DIR/shake256_vectors.hex"
set vec_dst "$PROJ_DIR/falcon_shake256.sim/sim_1/behav/xsim/shake256_vectors.hex"
file mkdir [file dirname $vec_dst]
file copy -force "$SIM/shake256_vectors.hex" $vec_dst

add_files -fileset sim_1 -norecurse "$SCRIPT_DIR/shake256_vectors.hex"
set_property file_type {Memory Initialization Files} \
    [get_files "$SCRIPT_DIR/shake256_vectors.hex"]

# ── Simulation settings ───────────────────────────────────────────────────────
set_property top     shake256_tb [get_filesets sim_1]
set_property top_lib work        [get_filesets sim_1]

set_property -name {xsim.simulate.runtime}     -value {all}  -objects [get_filesets sim_1]
set_property -name {xsim.simulate.log_all_signals} -value {true} -objects [get_filesets sim_1]

update_compile_order -fileset sources_1
update_compile_order -fileset sim_1

puts ""
puts "=================================================="
puts " Project: $PROJ_DIR/falcon_shake256.xpr"
puts "=================================================="
puts " Open in Vivado GUI, then:"
puts "   Flow Navigator → Simulation → Run Behavioral Simulation"
puts ""
puts " Key signals to add to waveform:"
puts "   shake256_tb/dut/fsm"
puts "   shake256_tb/dut/state[63:0]"
puts "   shake256_tb/dut/byte_pos"
puts "   shake256_tb/dut/sq_pos"
puts "   shake256_tb/dut/busy"
puts "   shake256_tb/dut/ready"
puts "=================================================="
