// fp_i2flt_int32_s.sv
// Signed 32-bit integer → IEEE 754 double, 1-cycle pipeline.
// Replaces Vivado Floating Point v7.1 IP (Int32→Double, NonBlocking, 1-cycle).
// Conversion is always exact: int32 has 31 significant bits ≤ 52 mantissa bits.

`timescale 1ns/1ps

module fp_i2flt_int32_s (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [31:0] s_axis_a_tdata,
    output wire        m_axis_result_tvalid,
    output reg  [63:0] m_axis_result_tdata
);
    assign m_axis_result_tvalid = 1'b1;

    wire        sign    = s_axis_a_tdata[31];
    wire [31:0] abs_val = sign ? (~s_axis_a_tdata + 32'd1) : s_axis_a_tdata;

    // Find position of highest set bit (k) in abs_val [31:0]
    // k determines the biased exponent: exp_out = 1023 + k
    function automatic [5:0] msb_pos32;
        input [31:0] v;
        integer i;
        begin
            msb_pos32 = 6'd0;
            for (i = 0; i <= 31; i = i + 1)
                if (v[i]) msb_pos32 = 6'(i);
        end
    endfunction

    wire [5:0]  k        = msb_pos32(abs_val);
    wire [10:0] exp_out  = 11'd1023 + {5'd0, k};
    // Align hidden bit to position 52: shift left by (52-k)
    wire [51:0] mant_out = (abs_val << (6'd52 - k)) >> 6'd0;  // lower 52 bits after shift

    // Build mantissa: shift abs_val left so bit k aligns to bit 52, then take [51:0]
    wire [63:0] shifted  = {32'd0, abs_val} << (6'd52 - k);
    wire [51:0] mant52   = shifted[51:0];

    wire [63:0] result_comb = (abs_val == 32'd0) ? 64'd0
                                                  : {sign, exp_out, mant52};

    always_ff @(posedge aclk) m_axis_result_tdata <= result_comb;

endmodule
