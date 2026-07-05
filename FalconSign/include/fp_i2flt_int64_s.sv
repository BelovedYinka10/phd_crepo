// fp_i2flt_int64_s.sv
// Signed 64-bit integer → IEEE 754 double, 1-cycle pipeline.
// Replaces Vivado Floating Point v7.1 IP (Int64→Double, NonBlocking, 1-cycle).
// Used in berexp.sv. Values up to 2^63 are handled; rounding is round-to-zero
// for values with more than 52 significant bits (acceptable for berexp comparison).

`timescale 1ns/1ps

module fp_i2flt_int64_s (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [63:0] s_axis_a_tdata,
    output wire        m_axis_result_tvalid,
    output reg  [63:0] m_axis_result_tdata
);
    assign m_axis_result_tvalid = 1'b1;

    wire        sign    = s_axis_a_tdata[63];
    // Absolute value (two's complement negate if negative)
    wire [63:0] abs_val = sign ? (~s_axis_a_tdata + 64'd1) : s_axis_a_tdata;

    // Find position of highest set bit k in abs_val[63:0]
    function automatic [6:0] msb_pos64;
        input [63:0] v;
        integer i;
        begin
            msb_pos64 = 7'd0;
            for (i = 0; i <= 63; i = i + 1)
                if (v[i]) msb_pos64 = 7'(i);
        end
    endfunction

    wire [6:0]  k       = msb_pos64(abs_val);       // position of MSB
    wire [10:0] exp_out = 11'd1023 + {4'd0, k};     // biased exponent

    // Align: shift abs_val so bit k is at position 52
    // k <= 52: shift left by (52-k)
    // k >  52: shift right by (k-52), truncating low bits (round to zero)
    wire        need_rshift = (k > 7'd52);
    wire [6:0]  lsh = need_rshift ? 7'd0       : (7'd52 - k);
    wire [6:0]  rsh = need_rshift ? (k - 7'd52) : 7'd0;

    wire [115:0] shifted_l = {52'd0, abs_val} << lsh;
    wire [63:0]  shifted_r = abs_val >> rsh;

    wire [51:0] mant52 = need_rshift ? shifted_r[51:0] : shifted_l[51:0];

    wire [63:0] result_comb = (abs_val == 64'd0) ? 64'd0
                                                  : {sign, exp_out, mant52};

    always_ff @(posedge aclk) m_axis_result_tdata <= result_comb;

endmodule
