// floating_point_sub_2pip.sv
// IEEE 754 double-precision subtraction (A - B), 2-cycle pipeline.
// Replaces Vivado Floating Point v7.1 IP (Low Latency, NonBlocking, 2-cycle).
// Implemented as addition with B's sign flipped.

`timescale 1ns/1ps

module floating_point_sub_2pip (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [63:0] s_axis_a_tdata,
    input  wire        s_axis_b_tvalid,
    input  wire [63:0] s_axis_b_tdata,
    output wire        m_axis_result_tvalid,
    output wire [63:0] m_axis_result_tdata
);
    // A - B  =  A + (-B)  →  flip B's sign bit and feed into the adder
    wire [63:0] b_neg = {~s_axis_b_tdata[63], s_axis_b_tdata[62:0]};

    floating_point_add_2pip u_add (
        .aclk                 (aclk),
        .s_axis_a_tvalid      (s_axis_a_tvalid),
        .s_axis_a_tdata       (s_axis_a_tdata),
        .s_axis_b_tvalid      (s_axis_b_tvalid),
        .s_axis_b_tdata       (b_neg),
        .m_axis_result_tvalid (m_axis_result_tvalid),
        .m_axis_result_tdata  (m_axis_result_tdata)
    );

endmodule
