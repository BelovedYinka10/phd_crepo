# create_project.tcl — Vivado 2019+ GUI project for FalconSign SamplerZ
#
# Run once to create the project:
#   vivado -mode batch -source create_project.tcl
#
# Then open falcon_sampler/falcon_sampler.xpr and click:
#   Flow Navigator → Simulation → Run Simulation → Run Behavioral Simulation
#
# Vivado will launch XSim and open the waveform viewer automatically.
# Add any signal from the hierarchy panel to the waveform to inspect it.

set SCRIPT_DIR [file normalize [file dirname [info script]]]
set ROOT       [file normalize "$SCRIPT_DIR/.."]
set INC        "$ROOT/include"
set SAM        "$ROOT/sampler"
set SIM        "$ROOT/sampler_sim"
set PROJ_DIR   "$SCRIPT_DIR/falcon_sampler"

# ── Create project ─────────────────────────────────────────────────────────────
# Part: ZCU104 (Zynq UltraScale+). Change if targeting a different board.
create_project falcon_sampler $PROJ_DIR \
    -part xczu7ev-ffvc1156-2-e -force

set_property default_lib       work          [current_project]
set_property target_language   SystemVerilog [current_project]
set_property simulator_language Mixed        [current_project]

# ── Design sources (packages first) ───────────────────────────────────────────
add_files -norecurse [list \
    $INC/falconsoar_pkg.sv \
    $INC/sample_pkg.sv \
    $INC/time_counter.sv \
    $INC/fpr_mul_opt_2pipe.sv \
    $INC/fpr_mul_const_opt_2pipe.sv \
    $INC/floating_point_add_2pip.sv \
    $INC/floating_point_sub_2pip.sv \
    $INC/fp_sub_s.sv \
    $INC/fp_mult_s.sv \
    $INC/fp_flt2i_int32_s.sv \
    $INC/fp_i2flt_int32_s.sv \
    $INC/fp_i2flt_int32.sv \
    $INC/fp_flt2i_int64_s.sv \
    $INC/fp_i2flt_int64_s.sv \
    $INC/counter_ce.sv \
    $INC/rbsh.sv \
    $SAM/chacha20.sv \
    $SAM/berexp.sv \
    $SAM/fpr_cal.sv \
    $SAM/pre_samp.sv \
    $SAM/samp_loop.sv \
    $SAM/refill_control.sv \
    $SAM/samplerz.sv \
    $SAM/samplerz_top.sv \
]

# Package files: mark as global include so every compilation unit sees them
foreach pkg_file [list $INC/falconsoar_pkg.sv $INC/sample_pkg.sv] {
    set_property file_type      {SystemVerilog Header} [get_files $pkg_file]
    set_property is_global_include true                [get_files $pkg_file]
}

# Add include directory so `include search works for other files
set_property include_dirs [list $INC] [current_fileset]

# ── Testbench (simulation fileset only) ───────────────────────────────────────
add_files -fileset sim_1 -norecurse $SIM/samplerz_tb.sv
set_property include_dirs [list $INC] [get_filesets sim_1]

# ── memory.hex — make sure it exists, then add as simulation data file ─────────
set mem_src "$SIM/memory.hex"

# Generate memory.hex if it doesn't exist
if {![file exists $mem_src]} {
    puts "Generating memory.hex..."
    exec python3 $SIM/gen_mem.py $mem_src
}

# Copy to project directory (XSim looks here during simulation)
file copy -force $mem_src $SCRIPT_DIR/memory.hex

# Also pre-populate the XSim working directory
set mem_dst "$PROJ_DIR/falcon_sampler.sim/sim_1/behav/xsim/memory.hex"
file mkdir [file dirname $mem_dst]
file copy -force $mem_src $mem_dst

# Add as a simulation file so Vivado knows about it
add_files -fileset sim_1 -norecurse $SCRIPT_DIR/memory.hex
set_property file_type {Memory Initialization Files} \
    [get_files $SCRIPT_DIR/memory.hex]

puts "memory.hex ready."

# ── Simulation settings ────────────────────────────────────────────────────────
set_property top     samplerz_tb [get_filesets sim_1]
set_property top_lib work        [get_filesets sim_1]

# Run until $finish (not a fixed time limit)
set_property -name {xsim.simulate.runtime} \
    -value {all} \
    -objects [get_filesets sim_1]

# Log all signals so the waveform database is complete
set_property -name {xsim.simulate.log_all_signals} \
    -value {true} \
    -objects [get_filesets sim_1]

# ── Save ───────────────────────────────────────────────────────────────────────
update_compile_order -fileset sources_1
update_compile_order -fileset sim_1

puts ""
puts "=================================================="
puts " Project ready: $PROJ_DIR/falcon_sampler.xpr"
puts "=================================================="
puts " Open in Vivado, then:"
puts "   Flow Navigator → Simulation"
puts "   → Run Simulation → Run Behavioral Simulation"
puts ""
puts " In the waveform viewer, add signals to inspect:"
puts "   samplerz_tb/dut/i_samplerz/i_samp_loop/random_bytes"
puts "   samplerz_tb/dut/i_samplerz/i_samp_loop/int_z"
puts "   samplerz_tb/dut/i_samplerz/i_berexp/fpr_x"
puts "   samplerz_tb/dut/i_samplerz/i_berexp/done"
puts "   samplerz_tb/dut/task_itf/op_done"
puts "=================================================="
