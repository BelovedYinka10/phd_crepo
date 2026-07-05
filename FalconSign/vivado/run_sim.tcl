# run_sim.tcl — FalconSign SamplerZ simulation in Vivado XSim
#
# Option A — batch mode (no GUI):
#   cd FalconSign/vivado
#   vivado -mode batch -source run_sim.tcl
#
# Option B — from Vivado Tcl console (GUI open):
#   source path/to/run_sim.tcl
#
# What this does:
#   1. Compiles all RTL + testbench with xvlog
#   2. Elaborates the design with xelab
#   3. Runs simulation to $finish and prints results
#   4. (Optional) Opens waveform if GUI is running
#
# No Vivado IP generation required — uses our RTL replacement modules
# for all floating-point operations (same files as the Verilator flow).
#
# Prerequisites:
#   1. Generate memory.hex first:
#        cd FalconSign/sampler_sim && python3 gen_mem.py memory.hex
#   2. Copy (or symlink) memory.hex into FalconSign/vivado/:
#        cp ../sampler_sim/memory.hex .

# ── Paths ─────────────────────────────────────────────────────────────────────
set SCRIPT_DIR [file dirname [file normalize [info script]]]
set ROOT       [file normalize "$SCRIPT_DIR/.."]
set INC        "$ROOT/include"
set SAM        "$ROOT/sampler"
set SIM        "$ROOT/sampler_sim"

# ── All source files (packages first, then modules, testbench last) ───────────
set SV_FILES [list \
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
    $SIM/samplerz_tb.sv \
]

# ── Step 1: compile ───────────────────────────────────────────────────────────
puts "=== Compiling SystemVerilog sources ==="
set xvlog_cmd "xvlog --sv --incr -L work"
foreach f $SV_FILES {
    append xvlog_cmd " \"$f\""
}
exec {*}[split $xvlog_cmd] >@stdout 2>@stderr

# ── Step 2: elaborate ─────────────────────────────────────────────────────────
puts "\n=== Elaborating: top = samplerz_tb ==="
exec xelab -debug typical samplerz_tb -snapshot samplerz_tb_snap \
    -log elaborate.log >@stdout 2>@stderr

# ── Step 3: simulate ──────────────────────────────────────────────────────────
puts "\n=== Running simulation ==="
if {[catch {
    exec xsim samplerz_tb_snap --runall -log simulate.log >@stdout 2>@stderr
} err]} {
    puts "Simulation finished (exit via \$finish): $err"
}

# ── Step 4: open waveform (GUI mode only) ────────────────────────────────────
if {[info exists ::env(DISPLAY)] || [string equal [lindex [split [info hostname] .] 0] ""] == 0} {
    catch {
        open_wave_database samplerz_tb_snap.wdb
        puts "\nWaveform database: samplerz_tb_snap.wdb"
        puts "Use 'open_wave_database samplerz_tb_snap.wdb' in the Tcl console to view."
    }
}

puts "\n=== Done. Check simulate.log for full output. ==="
