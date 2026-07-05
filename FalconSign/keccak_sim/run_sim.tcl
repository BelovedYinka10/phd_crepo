# run_sim.tcl — SHAKE256 batch simulation in Vivado XSim
#
# Batch mode:
#   cd FalconSign/keccak_sim
#   python3 gen_vectors.py          # generate shake256_vectors.hex first
#   vivado -mode batch -source run_sim.tcl

set SCRIPT_DIR [file dirname [file normalize [info script]]]
set ROOT       [file normalize "$SCRIPT_DIR/.."]
set KECCAK     "$ROOT/keccak"
set SIM        "$SCRIPT_DIR"

set SV_FILES [list \
    $KECCAK/keccak_f1600.sv \
    $KECCAK/shake256.sv     \
    $SIM/shake256_tb.sv     \
]

# ── Step 1: compile ───────────────────────────────────────────────────────────
puts "=== Compiling ==="
exec xvlog --sv --nolog {*}$SV_FILES >@stdout 2>@stderr

# ── Step 2: elaborate ─────────────────────────────────────────────────────────
puts "\n=== Elaborating: top = shake256_tb ==="
exec xelab -debug typical shake256_tb -snapshot shake256_tb_snap \
    --nolog >@stdout 2>@stderr

# ── Step 3: simulate ──────────────────────────────────────────────────────────
puts "\n=== Simulating ==="
if {[catch {
    exec xsim shake256_tb_snap --runall --nolog >@stdout 2>@stderr
} err]} {
    puts "Simulation ended: $err"
}

puts "\n=== Done ==="
