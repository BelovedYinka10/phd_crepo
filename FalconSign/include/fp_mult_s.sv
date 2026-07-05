// fp_mult_s.sv
// IEEE 754 double-precision multiplication, 2-cycle pipeline.
// Replaces Vivado Floating Point v7.1 IP (High Speed, NonBlocking, 2-cycle).
// Used in the SamplerZ datapath (fpr_cal, berexp).
// Wraps fpr_mul_opt_2pipe with the full AXI4-S tvalid port set.

`timescale 1ns/1ps

module fp_mult_s (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [63:0] s_axis_a_tdata,
    input  wire        s_axis_b_tvalid,
    input  wire [63:0] s_axis_b_tdata,
    output wire        m_axis_result_tvalid,
    output wire [63:0] m_axis_result_tdata
);
    assign m_axis_result_tvalid = 1'b1;

    fpr_mul_opt_2pipe u_mul (
        .aclk                (aclk),
        .s_axis_a_tdata      (s_axis_a_tdata),
        .s_axis_b_tdata      (s_axis_b_tdata),
        .m_axis_result_tdata (m_axis_result_tdata)
    );

endmodule
