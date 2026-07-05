# run_vivado.tcl — batch synth + impl of the BaseSampler, report Fmax + area.
#
# Usage (from this directory):
#   vivado -mode batch -source run_vivado.tcl
#
# Edit PART to your target device before running.  Outputs:
#   util.rpt     — LUT/FF/etc. utilization
#   timing.rpt   — full timing summary
#   and Fmax printed to the console + fmax.rpt.

# ---- target part: CHANGE THIS to your board/device ----------------------
set PART   xc7a100tcsg324-1   ;# Artix-7 example; e.g. Zynq xc7z020clg400-1,
                              ;#  UltraScale+ xczu3eg-sbva484-1-e, etc.
set TOP    base_sampler_top
set PERIOD 2.000              ;# must match create_clock in base_sampler.xdc (ns)

# ---- read sources -------------------------------------------------------
read_verilog ../base_sampler.v
read_verilog base_sampler_top.v
read_xdc     base_sampler.xdc

# ---- synthesis ----------------------------------------------------------
synth_design -top $TOP -part $PART
opt_design

# ---- implementation -----------------------------------------------------
place_design
route_design

# ---- reports ------------------------------------------------------------
report_utilization      -file util.rpt
report_timing_summary   -file timing.rpt -delay_type max -max_paths 10

# ---- compute Fmax from worst setup slack --------------------------------
set wns [get_property SLACK [get_timing_paths -setup -max_paths 1 -nworst 1]]
if {$wns eq ""} { set wns 0.0 }
set achieved [expr {$PERIOD - $wns}]
set fmax     [expr {1000.0 / $achieved}]

set fh [open fmax.rpt w]
puts $fh "part        = $PART"
puts $fh "clk target  = $PERIOD ns ([format %.1f [expr {1000.0/$PERIOD}]] MHz)"
puts $fh "WNS         = $wns ns"
puts $fh "achieved    = [format %.3f $achieved] ns"
puts $fh "Fmax        = [format %.1f $fmax] MHz"
close $fh

puts "=================================================="
puts "  BaseSampler datapath (compare + popcount)"
puts "  part   : $PART"
puts "  WNS    : $wns ns  (target period $PERIOD ns)"
puts "  Fmax   : [format %.1f $fmax] MHz   => 1 sample / cycle"
puts "  (PRNG feed not included; it will lower the system Fmax)"
puts "  see util.rpt for LUT/FF area"
puts "=================================================="
