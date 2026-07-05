// rbsh.sv
// Right barrel (circular) shift — combinational.
// Interface matches DesignWare DW_rbsh.
//
//   SH_TC = 0 : SH is unsigned → right circular shift A by SH positions
//   SH_TC = 1 : SH is signed   → positive SH: right shift; negative SH: left shift
//
// All callers in this project use SH_TC = 0.
//
// Right circular shift by k : B[i] = A[(i+k) % A_width]
//   Implemented via the doubling trick: {A,A}[A_width+k-1 : k]
//
// Variable part-select [msb -: width] is supported by Vivado / ModelSim.

`timescale 1ns/1ps

module rbsh #(
    parameter int A_width  = 8,
    parameter int SH_width = 3
) (
    input  wire [A_width-1:0]  A,
    input  wire [SH_width-1:0] SH,
    input  wire                SH_TC,  // 1 = signed shift amount
    output wire [A_width-1:0]  B
);

    // Effective unsigned right-shift amount (mod A_width)
    wire                  sh_neg = SH_TC & SH[SH_width-1];
    wire [SH_width-1:0]   sh_abs = sh_neg ? (-SH) : SH;

    // When SH_TC=1 and SH is negative, a left circular shift of sh_abs
    // equals a right circular shift of (A_width - sh_abs).
    // We reduce both to a single right-rotation amount in [0, A_width).
    localparam int AMSB = $clog2(A_width) + 1;  // enough bits for A_width value
    wire [AMSB-1:0] sh_right = sh_neg
                               ? (AMSB'(A_width) - AMSB'(sh_abs))
                               : AMSB'(sh_abs);

    // Modulo A_width (only needed if caller can exceed A_width; in practice
    // rotate_num = stage % log2(n) which is always < A_width)
    wire [AMSB-1:0] sh_eff = (sh_right >= AMSB'(A_width)) ? (sh_right - AMSB'(A_width))
                                                            :  sh_right;

    // Double A so we can extract any rotation as a contiguous slice
    wire [2*A_width-1:0] doubled = {A, A};

    // B = doubled[A_width + sh_eff - 1 : sh_eff]  (descending part-select)
    assign B = doubled[A_width + sh_eff - 1 -: A_width];

endmodule