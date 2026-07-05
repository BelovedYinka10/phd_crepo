// fp_i2flt_int32.sv
// Signed 32-bit integer → IEEE 754 double, 1-cycle pipeline.
// Used in samplerz.sv to convert the final sample value to double.
// Differs from fp_i2flt_int32_s only in the port list:
//   - has m_axis_result_tready (input, ignored) instead of m_axis_result_tvalid

`timescale 1ns/1ps

module fp_i2flt_int32 (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [31:0] s_axis_a_tdata,
    input  wire        m_axis_result_tready,   // ignored — always ready
    output reg  [63:0] m_axis_result_tdata
);
    wire        sign    = s_axis_a_tdata[31];
    wire [31:0] abs_val = sign ? (~s_axis_a_tdata + 32'd1) : s_axis_a_tdata;

    function automatic [5:0] msb_pos32;
        input [31:0] v;
        integer i;
        begin
            msb_pos32 = 6'd0;
            for (i = 0; i <= 31; i = i + 1)
                if (v[i]) msb_pos32 = 6'(i);
        end
    endfunction

    wire [5:0]  k       = msb_pos32(abs_val);
    wire [10:0] exp_out = 11'd1023 + {5'd0, k};
    wire [63:0] shifted = {32'd0, abs_val} << (6'd52 - k);
    wire [51:0] mant52  = shifted[51:0];

    wire [63:0] result_comb = (abs_val == 32'd0) ? 64'd0
                                                  : {sign, exp_out, mant52};

    always_ff @(posedge aclk) m_axis_result_tdata <= result_comb;

endmodule
