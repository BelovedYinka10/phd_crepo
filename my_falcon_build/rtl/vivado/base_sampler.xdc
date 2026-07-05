# base_sampler.xdc — timing constraint for Fmax exploration.
#
# Set an aggressive target period, then read the achieved Fmax from the
# implemented design's WNS:  Fmax = 1000 / (PERIOD_ns - WNS_ns)  [MHz].
#
# Tighten PERIOD until WNS goes slightly negative to find the true ceiling,
# or just read Fmax off one run (the run_vivado.tcl script computes it).
#
# 2.000 ns => 500 MHz target.  Adjust to taste / to your part's capability.
create_clock -name clk -period 2.000 [get_ports clk]

# This is a pure datapath block with no real I/O timing budget to model;
# relax I/O delays so the report reflects the internal logic path, not pad
# timing.  (For an OOC / sub-block characterisation this is the honest view.)
set_input_delay  -clock clk 0.000 [get_ports {u_in[*] rst}]
set_output_delay -clock clk 0.000 [get_ports {z0_out[*]}]
