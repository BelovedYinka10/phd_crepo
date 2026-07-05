// fp_flt2i_int32_s.sv
// IEEE 754 double → signed 32-bit integer, 1-cycle pipeline.
// Replaces Vivado Floating Point v7.1 IP (Double→Int32, NonBlocking, 1-cycle).
// Rounding mode: truncate toward zero (matches C cast behaviour).
// Saturates on overflow: >2^31-1 → INT_MAX, <-2^31 → INT_MIN.

`timescale 1ns/1ps

module fp_flt2i_int32_s (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [63:0] s_axis_a_tdata,
    output wire        m_axis_result_tvalid,
    output reg  [31:0] m_axis_result_tdata
);
    assign m_axis_result_tvalid = 1'b1;

    wire        sign = s_axis_a_tdata[63];
    wire [10:0] exp  = s_axis_a_tdata[62:52];
    wire [51:0] mant = s_axis_a_tdata[51:0];

    wire is_nan  = (exp == 11'h7FF) && (mant != 52'd0);
    wire is_inf  = (exp == 11'h7FF) && (mant == 52'd0);
    wire is_zero = (exp == 11'd0);

    // Actual exponent k = exp - 1023
    // k < 0  (exp < 1023) : |value| < 1 → truncate to 0
    // k >= 31 (exp >= 1054): overflow for int32
    wire        k_neg  = (exp < 11'd1023);
    wire        k_ovf  = (exp >= 11'd1054);
    wire [10:0] k_val  = exp - 11'd1023;   // actual exponent; valid when !k_neg && !k_ovf

    // k in [0..30]: shift significand right by (52-k) to extract integer
    wire [5:0]  rshift    = 6'd52 - k_val[5:0];  // k<=30, so k_val[5:0] is exact
    wire [52:0] sig       = {1'b1, mant};
    wire [52:0] abs_full  = sig >> rshift;
    wire [30:0] abs_val   = abs_full[30:0];

    wire [31:0] pos_result = {1'b0, abs_val};
    wire [31:0] neg_result = ~{1'b0, abs_val} + 32'd1;
    wire [31:0] normal     = sign ? neg_result : pos_result;

    wire [31:0] result_comb =
        (is_nan | is_inf | k_ovf) ? (sign ? 32'h8000_0000 : 32'h7FFF_FFFF) :
        (is_zero | k_neg)         ? 32'd0 :
                                    normal;

    always_ff @(posedge aclk) m_axis_result_tdata <= result_comb;

endmodule
