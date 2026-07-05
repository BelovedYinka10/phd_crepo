// fpr_mul_const_opt_2pipe.sv
// IEEE 754 double-precision multiply by 0.5 (FPR_HALF = 64'h3FE0000000000000).
// 2-cycle registered pipeline to match the latency of fpr_mul_opt_2pipe.
//
// Interface mirrors the Xilinx Floating Point IP AXI4-S wrapper used throughout
// bfu_fft.sv (aclk / s_axis_a_tdata / m_axis_result_tdata).
//
// Multiply by 0.5 reduces to decrementing the biased exponent by 1:
//
//   normal (exp > 1) : {sign, exp-1, mant}          -- exact, no rounding needed
//   exp == 1         : {sign, 0,     {1,mant}>>1}    -- smallest normal → denormal
//   zero / denormal  : unchanged (result is 0 or tiny denormal / 2)
//   inf              : unchanged
//   NaN              : unchanged (quiet NaN propagated)
//
// In Falcon, all FFT operands are well within normal range, so the exp==1 and
// denormal branches are never exercised in practice.

`timescale 1ns/1ps

module fpr_mul_const_opt_2pipe (
    input  wire         aclk,
    input  wire  [63:0] s_axis_a_tdata,
    output reg   [63:0] m_axis_result_tdata
);

    // IEEE 754 double fields
    wire        sign_a = s_axis_a_tdata[63];
    wire [10:0] exp_a  = s_axis_a_tdata[62:52];
    wire [51:0] mant_a = s_axis_a_tdata[51:0];

    wire is_nan  = (exp_a == 11'h7FF) && (mant_a != 52'd0);
    wire is_inf  = (exp_a == 11'h7FF) && (mant_a == 52'd0);
    wire is_zero = (exp_a == 11'd0)   && (mant_a == 52'd0);

    // exp == 1: normal → denormal after ×0.5
    // Shift the full significand (hidden bit 1 concatenated) right by 1.
    wire [51:0] denorm_mant = {1'b1, mant_a[51:1]};  // {1, mant} >> 1 (drop LSB)

    // Combinational result
    wire [63:0] result_comb =
        (is_nan | is_inf | is_zero)  ? s_axis_a_tdata :                // pass-through
        (exp_a == 11'd1)             ? {sign_a, 11'd0, denorm_mant} :  // → denormal
                                       {sign_a, exp_a - 11'd1, mant_a};// normal ×0.5

    // Two-stage pipeline
    reg [63:0] pipe1;
    always_ff @(posedge aclk) pipe1               <= result_comb;
    always_ff @(posedge aclk) m_axis_result_tdata <= pipe1;

endmodule