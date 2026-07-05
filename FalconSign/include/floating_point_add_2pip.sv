// floating_point_add_2pip.sv
// IEEE 754 double-precision addition, 2-cycle pipeline.
// Replaces Vivado Floating Point v7.1 IP (Low Latency, NonBlocking, 2-cycle).
// tvalid ports accepted but ignored — result is always available after 2 cycles.
//
// Stage 1: extract fields, detect specials, align significands
// Stage 2: add significands, normalize, round, pack

`timescale 1ns/1ps

module floating_point_add_2pip (
    input  wire        aclk,
    input  wire        s_axis_a_tvalid,
    input  wire [63:0] s_axis_a_tdata,
    input  wire        s_axis_b_tvalid,
    input  wire [63:0] s_axis_b_tdata,
    output wire        m_axis_result_tvalid,
    output reg  [63:0] m_axis_result_tdata
);
    assign m_axis_result_tvalid = 1'b1;

    // ── Field extraction ────────────────────────────────────────────────────────
    wire        sa = s_axis_a_tdata[63];
    wire [10:0] ea = s_axis_a_tdata[62:52];
    wire [51:0] ma = s_axis_a_tdata[51:0];
    wire        sb = s_axis_b_tdata[63];
    wire [10:0] eb = s_axis_b_tdata[62:52];
    wire [51:0] mb = s_axis_b_tdata[51:0];

    wire [52:0] sig_a = {ea != 11'd0, ma};
    wire [52:0] sig_b = {eb != 11'd0, mb};

    wire nan_a  = (ea == 11'h7FF) && (ma != 52'd0);
    wire nan_b  = (eb == 11'h7FF) && (mb != 52'd0);
    wire inf_a  = (ea == 11'h7FF) && (ma == 52'd0);
    wire inf_b  = (eb == 11'h7FF) && (mb == 52'd0);
    wire zero_a = (ea == 11'd0);
    wire zero_b = (eb == 11'd0);

    // ── Alignment ────────────────────────────────────────────────────────────────
    wire a_larger   = (ea > eb) | ((ea == eb) & (sig_a >= sig_b));
    wire [10:0] exp_l  = a_larger ? ea    : eb;
    wire [52:0] sig_l  = a_larger ? sig_a : sig_b;
    wire        sign_l = a_larger ? sa    : sb;
    wire [10:0] exp_s  = a_larger ? eb    : ea;
    wire [52:0] sig_s  = a_larger ? sig_b : sig_a;
    wire        sign_s = a_larger ? sb    : sa;

    wire [6:0] exp_diff = (exp_l >= exp_s) ? (7'(exp_l) - 7'(exp_s)) : 7'd0;
    wire [6:0] shift    = (exp_diff > 7'd54) ? 7'd54 : exp_diff;
    wire [52:0] sig_s_sh = sig_s >> shift;

    wire eff_sub    = (sign_l != sign_s);
    wire nan_out    = nan_a | nan_b | (inf_a & inf_b & eff_sub);
    wire inf_out    = (inf_a | inf_b) & ~nan_out;
    wire inf_sign   = inf_a ? sa : sb;
    wire zero_out   = zero_a & zero_b;

    // ── Stage 1 registers ─────────────────────────────────────────────────────
    reg         p1_sign, p1_eff_sub, p1_nan, p1_inf, p1_zero, p1_inf_sign;
    reg [11:0]  p1_exp;
    reg [52:0]  p1_sig_l, p1_sig_s_sh;

    always_ff @(posedge aclk) begin
        p1_sign     <= sign_l;
        p1_eff_sub  <= eff_sub;
        p1_exp      <= {1'b0, exp_l};
        p1_sig_l    <= sig_l;
        p1_sig_s_sh <= sig_s_sh;
        p1_nan      <= nan_out;
        p1_inf      <= inf_out;
        p1_zero     <= zero_out;
        p1_inf_sign <= inf_sign;
    end

    // ── Stage 2: add/sub, normalize, pack ───────────────────────────────────────
    // Leading-zero count for 53-bit significand (for normalization after cancellation)
    function automatic [5:0] lzc53;
        input [52:0] v;
        integer i;
        begin
            lzc53 = 6'd53;  // all-zero case
            for (i = 0; i <= 52; i = i + 1)
                if (v[i]) lzc53 = 6'(52 - i);  // iterate LSB→MSB; last hit = MSB position
        end
    endfunction

    wire [53:0] sum_raw = p1_eff_sub ? ({1'b0,p1_sig_l} - {1'b0,p1_sig_s_sh})
                                     : ({1'b0,p1_sig_l} + {1'b0,p1_sig_s_sh});
    wire        carry   = sum_raw[53];
    wire [5:0]  lz      = lzc53(sum_raw[52:0]);
    wire [52:0] sig_shifted = sum_raw[52:0] << lz;
    wire [51:0] mant_out = carry ? sum_raw[52:1] : sig_shifted[51:0];
    wire [11:0] exp_out  = carry            ? (p1_exp + 12'd1) :
                           (sum_raw == '0)  ? 12'd0 :
                                              (p1_exp - {6'd0, lz});
    wire ovf  = ~exp_out[11] && (exp_out >= 12'd2047);
    wire undf = exp_out[11] || (exp_out == 12'd0 && sum_raw != '0);

    localparam [63:0] QNAN = 64'h7FF8_0000_0000_0000;

    always_ff @(posedge aclk) begin
        if      (p1_nan)                             m_axis_result_tdata <= QNAN;
        else if (p1_inf)                             m_axis_result_tdata <= {p1_inf_sign, 11'h7FF, 52'd0};
        else if (p1_zero || sum_raw == '0)           m_axis_result_tdata <= {p1_sign, 63'd0};
        else if (ovf)                                m_axis_result_tdata <= {p1_sign, 11'h7FF, 52'd0};
        else if (undf)                               m_axis_result_tdata <= {p1_sign, 63'd0};
        else                                         m_axis_result_tdata <= {p1_sign, exp_out[10:0], mant_out};
    end

endmodule