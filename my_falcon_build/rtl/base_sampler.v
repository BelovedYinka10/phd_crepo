`timescale 1ns/1ps
// base_sampler.v
//
// Falcon BaseSampler — half-Gaussian (sigma_max = 1.8205) via RCDT inversion.
// ==========================================================================
//
//   z0 = #{ k : u < RCDT[k] },   k = 0 .. 17
//
// The software reference (my_falcon_build/samplerz.py, basesampler()) must
// scan all 18 table entries sequentially to stay constant-time.  In hardware
// the same "touch everything, no early-exit" property is free: the 18
// comparisons happen in parallel in a single combinational cone, and a
// popcount adder tree reduces them to the sample.
//
//   72-bit u ──┬─[u<RCDT0 ]─┐
//              ├─[u<RCDT1 ]─┤
//              │     ...     ├─ popcount ─► z0 (0..18, 5 bits)
//              └─[u<RCDT17]─┘
//
// One sample per evaluation.  For N-lane throughput (the SIMD paper's
// "vectorization"), instantiate N copies sharing a wider PRNG.
//
// Bit-for-bit identical to samplerz.py basesampler():
//   u  = int.from_bytes(randombytes(9), "little")   // 72-bit LE word
//   z0 = sum(int(u < elt) for elt in RCDT)
//
module base_sampler (
    input  wire [71:0] u,    // 72-bit uniform random word
    output wire [4:0]  z0    // sample in 0 .. 18  (needs 5 bits)
);

    // RCDT[k] = round( Pr[X > k] * 2^72 )  for  X ~ half-Gaussian(sigma_max).
    // Values copied verbatim from samplerz.py RCDT[] (decimal).
    localparam [71:0]
        RCDT0  = 72'd3024686241123004913666,
        RCDT1  = 72'd1564742784480091954050,
        RCDT2  = 72'd636254429462080897535,
        RCDT3  = 72'd199560484645026482916,
        RCDT4  = 72'd47667343854657281903,
        RCDT5  = 72'd8595902006365044063,
        RCDT6  = 72'd1163297957344668388,
        RCDT7  = 72'd117656387352093658,
        RCDT8  = 72'd8867391802663976,
        RCDT9  = 72'd496969357462633,
        RCDT10 = 72'd20680885154299,
        RCDT11 = 72'd638331848991,
        RCDT12 = 72'd14602316184,
        RCDT13 = 72'd247426747,
        RCDT14 = 72'd3104126,
        RCDT15 = 72'd28824,
        RCDT16 = 72'd198,
        RCDT17 = 72'd1;

    // gt[k] = 1  iff  u < RCDT[k]   (unsigned 72-bit compare)
    wire [17:0] gt;
    assign gt[0]  = (u < RCDT0 );
    assign gt[1]  = (u < RCDT1 );
    assign gt[2]  = (u < RCDT2 );
    assign gt[3]  = (u < RCDT3 );
    assign gt[4]  = (u < RCDT4 );
    assign gt[5]  = (u < RCDT5 );
    assign gt[6]  = (u < RCDT6 );
    assign gt[7]  = (u < RCDT7 );
    assign gt[8]  = (u < RCDT8 );
    assign gt[9]  = (u < RCDT9 );
    assign gt[10] = (u < RCDT10);
    assign gt[11] = (u < RCDT11);
    assign gt[12] = (u < RCDT12);
    assign gt[13] = (u < RCDT13);
    assign gt[14] = (u < RCDT14);
    assign gt[15] = (u < RCDT15);
    assign gt[16] = (u < RCDT16);
    assign gt[17] = (u < RCDT17);

    // popcount(gt) -> z0  (synthesises to an adder tree)
    function [4:0] popcount18;
        input [17:0] b;
        integer i;
        reg [4:0] acc;
        begin
            acc = 5'd0;
            for (i = 0; i < 18; i = i + 1)
                acc = acc + b[i];
            popcount18 = acc;
        end
    endfunction

    assign z0 = popcount18(gt);

endmodule
