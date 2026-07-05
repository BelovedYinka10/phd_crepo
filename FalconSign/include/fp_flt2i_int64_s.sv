// fp_flt2i_int64_s.sv
// IEEE 754 double → signed 64-bit integer, 1-cycle pipeline.
// Replaces Vivado Floating Point v7.1 IP (Double→Int64, NonBlocking, 1-cycle).
// Used in berexp.sv for the approxexp bit-comparison against 2^63·ccs·exp(-x).
// Rounding: truncate toward zero. Saturates on overflow.

`timescale 1ns/1ps

module fp_flt2i_int64_s (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [63:0] s_axis_a_tdata,
    output wire        m_axis_result_tvalid,
    output reg  [63:0] m_axis_result_tdata
);
    assign m_axis_result_tvalid = 1'b1;

    wire        sign = s_axis_a_tdata[63];
    wire [10:0] exp  = s_axis_a_tdata[62:52];
    wire [51:0] mant = s_axis_a_tdata[51:0];

    wire is_nan  = (exp == 11'h7FF) && (mant != 52'd0);
    wire is_inf  = (exp == 11'h7FF) && (mant == 52'd0);
    wire is_zero = (exp == 11'd0);

    // k = exp - 1023 (actual exponent)
    wire        k_neg = (exp < 11'd1023);          // |value| < 1 → 0
    wire        k_ovf = (exp >= 11'd1086);         // k >= 63 → overflow for int64

    wire [6:0]  k      = 7'(exp) - 7'd1023;       // k in [0..62] when valid
    wire [52:0] sig    = {1'b1, mant};

    // If k <= 52: right-shift sig by (52-k) to get integer
    // If k >  52: left-shift sig by (k-52) to get integer
    wire        lshift_mode = (k > 7'd52) && !k_neg && !k_ovf;
    wire [6:0]  rshift_amt  = 7'd52 - k;           // valid when k <= 52
    wire [6:0]  lshift_amt  = k - 7'd52;           // valid when k > 52

    wire [115:0] sig_shifted_l = {63'd0, sig} << lshift_amt;
    wire [52:0]  sig_shifted_r = sig >> rshift_amt;

    wire [63:0] abs_val = lshift_mode ? sig_shifted_l[63:0]
                                      : {11'd0, sig_shifted_r};

    wire [63:0] pos_result = abs_val;
    wire [63:0] neg_result = ~abs_val + 64'd1;
    wire [63:0] normal     = sign ? neg_result : pos_result;

    localparam [63:0] INT64_MAX = 64'h7FFF_FFFF_FFFF_FFFF;
    localparam [63:0] INT64_MIN = 64'h8000_0000_0000_0000;

    wire [63:0] result_comb =
        (is_nan | is_inf | k_ovf) ? (sign ? INT64_MIN : INT64_MAX) :
        (is_zero | k_neg)         ? 64'd0 :
                                    normal;

    always_ff @(posedge aclk) m_axis_result_tdata <= result_comb;

endmodule
