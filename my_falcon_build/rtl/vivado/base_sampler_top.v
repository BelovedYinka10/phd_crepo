`timescale 1ns/1ps
// base_sampler_top.v
//
// Registered wrapper around the combinational base_sampler core, so Vivado
// has a clocked path to measure (Fmax = max frequency of the input-FF ->
// comparators -> popcount -> output-FF path).
//
// Synthesise/implement THIS as top (not base_sampler.v alone), with a clock
// constraint on `clk` (see base_sampler.xdc).  The registers isolate the
// compare+popcount logic as the critical path.
//
// NOTE: this measures ONLY the comparator/popcount block.  The real sampler
// also needs a PRNG (SHAKE/ChaCha) to produce `u` — that block, not this one,
// sets the true system Fmax and dominates area.  Treat this Fmax as an upper
// bound on the BaseSampler datapath.
//
module base_sampler_top (
    input  wire        clk,
    input  wire        rst,
    input  wire [71:0] u_in,    // 72-bit uniform random word
    output reg  [4:0]  z0_out   // sampled z0 (0..18)
);

    reg  [71:0] u_q;
    wire [4:0]  z0_w;

    // Input register
    always @(posedge clk) begin
        if (rst) u_q <= 72'd0;
        else     u_q <= u_in;
    end

    // Combinational core under test
    base_sampler core (
        .u  (u_q),
        .z0 (z0_w)
    );

    // Output register
    always @(posedge clk) begin
        if (rst) z0_out <= 5'd0;
        else     z0_out <= z0_w;
    end

endmodule
