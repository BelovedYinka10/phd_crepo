// keccak_shake256_wrapper.sv — byte-stream SHAKE256 adapter around cshake-core's keccak_top.
//
// keccak_top (cshake-core/verilog/keccak_top.v, Yale/Jungk, validated standalone via
// cshake-core/Vivado) speaks a 32-bit-word, block-oriented protocol, not a byte stream:
//   1 command word   {cshake=0, mux256=1, out_len_bits[29:0]}     (once per hash)
//   per 136-byte block:
//     1 length word   {eof, 30'b0, block_len_bits[10:0]}
//     34 data words   (136 bytes; zero-padded past the real message on the last block)
//   then a fixed-length squeeze of out_len_bits/8 bytes.
// Both total lengths must be known before absorption starts, so this wrapper adds
// msg_len_bytes / out_len_bytes ports (latched on `init`) that shake256.sv did not need.
// Everything else (absorb_valid/absorb_data/absorb_last/ready, squeeze_ready/
// squeeze_data/squeeze_valid, busy) matches shake256.sv's port names so callers only
// need to add the two length ports.
//
// Byte<->word packing matches the cshake-core test-vector generator convention
// (testvectors_generator/generator.cpp): byte i -> word bits [8*(i%4) +: 8], i.e. the
// first byte of a 4-byte group lands in the word's low byte.
//
// out_len_bytes must be a multiple of 4 (round up to the nearest word — hash_to_point
// picks n*3 already, well clear of what it will ever actually consume).

`timescale 1ns/1ps

module keccak_shake256_wrapper (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        init,             // 1-cycle pulse: begin new hash
    input  wire [31:0] msg_len_bytes,    // total bytes to be absorbed (salt||message), valid at init
    input  wire [31:0] out_len_bytes,    // total bytes to be squeezed, valid at init (must be mult of 4)
    // absorb (byte stream)
    input  wire        absorb_valid,
    input  wire [7:0]  absorb_data,
    input  wire        absorb_last,      // sanity pulse only; length is already known via msg_len_bytes
    output wire        ready,            // 1 = accepting a message byte this cycle
    // squeeze (byte stream)
    input  wire        squeeze_ready,
    output wire [7:0]  squeeze_data,
    output wire        squeeze_valid,
    output wire        busy
);

    localparam RATE_BYTES = 136;   // SHAKE256 rate: 1088 bits

    // ── keccak_top instance ──────────────────────────────────────────────────
    wire        kt_rst;
    reg         kt_din_valid;
    wire        kt_din_ready;
    reg  [31:0] kt_din;
    wire        kt_dout_valid;
    reg         kt_dout_ready;
    wire [31:0] kt_dout;

    keccak_top u_keccak_top (
        .clk        (clk),
        .rst        (kt_rst),
        .din_valid  (kt_din_valid),
        .din_ready  (kt_din_ready),
        .din        (kt_din),
        .dout_valid (kt_dout_valid),
        .dout_ready (kt_dout_ready),
        .dout       (kt_dout)
    );

    // ── Wrapper FSM ───────────────────────────────────────────────────────────
    localparam W_IDLE     = 3'd0;
    localparam W_CMD      = 3'd1;
    localparam W_LEN      = 3'd2;
    localparam W_FILL     = 3'd3;
    localparam W_SENDWORD = 3'd4;
    localparam W_SQUEEZE  = 3'd5;

    reg [2:0]  wfsm;
    reg [31:0] total_len, total_out;
    reg [31:0] bytes_left;       // bytes not yet assigned to a block
    reg [31:0] msg_consumed;     // running count of real bytes packed so far
    reg [31:0] out_produced;     // running count of bytes delivered to squeeze_data
    reg [5:0]  word_idx;         // 0..33 within current block
    reg [1:0]  byte_idx;         // 0..3 within word_buf
    reg [31:0] word_buf;
    reg        block_eof_reg;
    reg [31:0] block_bytes_reg;

    wire [31:0] block_bytes_c = (bytes_left > RATE_BYTES) ? 32'd136 : bytes_left;
    wire        block_eof_c   = (bytes_left <= RATE_BYTES);

    wire [31:0] bytes_left_after   = bytes_left - block_bytes_reg;
    wire [31:0] next_block_bytes_c = (bytes_left_after > RATE_BYTES) ? 32'd136 : bytes_left_after;
    wire        next_block_eof_c   = (bytes_left_after <= RATE_BYTES);

    assign kt_rst = !rst_n;
    assign busy   = (wfsm != W_IDLE);
    assign ready  = (wfsm == W_FILL) && (msg_consumed < total_len);

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            wfsm          <= W_IDLE;
            kt_din_valid  <= 0;
            word_idx      <= 0;
            byte_idx      <= 0;
            msg_consumed  <= 0;
        end else begin
            case (wfsm)

            W_IDLE: begin
                kt_din_valid <= 0;
                if (init) begin
                    total_len    <= msg_len_bytes;
                    total_out    <= out_len_bytes;
                    bytes_left   <= msg_len_bytes;
                    msg_consumed <= 0;
                    kt_din       <= {1'b0, 1'b1, out_len_bytes[26:0], 3'b000}; // cshake=0, mux256=1, out_len_bytes*8
                    kt_din_valid <= 1;
                    wfsm         <= W_CMD;
                end
            end

            W_CMD: begin
                if (kt_din_valid && kt_din_ready) begin
                    kt_din_valid    <= 1;
                    kt_din          <= {block_eof_c, 20'b0, block_bytes_c[7:0], 3'b000};
                    wfsm            <= W_LEN;
                end
            end

            W_LEN: begin
                if (kt_din_valid && kt_din_ready) begin
                    block_eof_reg   <= block_eof_c;
                    block_bytes_reg <= block_bytes_c;
                    kt_din_valid    <= 0;
                    word_idx        <= 0;
                    byte_idx        <= 0;
                    wfsm            <= W_FILL;
                end
            end

            W_FILL: begin
                if (msg_consumed < total_len) begin
                    if (absorb_valid) begin
                        if (byte_idx == 2'd3) begin
                            kt_din       <= {absorb_data, word_buf[23:0]};
                            kt_din_valid <= 1;
                            wfsm         <= W_SENDWORD;
                        end else begin
                            word_buf[byte_idx*8 +: 8] <= absorb_data;
                            byte_idx <= byte_idx + 2'd1;
                        end
                        msg_consumed <= msg_consumed + 1'b1;
                    end
                end else begin
                    if (byte_idx == 2'd3) begin
                        kt_din       <= {8'h00, word_buf[23:0]};
                        kt_din_valid <= 1;
                        wfsm         <= W_SENDWORD;
                    end else begin
                        word_buf[byte_idx*8 +: 8] <= 8'h00;
                        byte_idx <= byte_idx + 2'd1;
                    end
                end
            end

            W_SENDWORD: begin
                if (kt_din_valid && kt_din_ready) begin
                    kt_din_valid <= 0;
                    byte_idx     <= 0;
                    if (word_idx == 6'd33) begin
                        if (block_eof_reg) begin
                            wfsm <= W_SQUEEZE;
                            out_produced <= 0;
                        end else begin
                            bytes_left   <= bytes_left_after;
                            kt_din_valid <= 1;
                            kt_din       <= {next_block_eof_c, 20'b0, next_block_bytes_c[7:0], 3'b000};
                            wfsm         <= W_LEN;
                        end
                    end else begin
                        word_idx <= word_idx + 6'd1;
                        wfsm     <= W_FILL;
                    end
                end
            end

            W_SQUEEZE: begin
                // byte unpacking/draining happens in the squeeze always_ff below;
                // once every requested byte has been delivered, go idle so the next
                // `init` can start a new hash.
                if (out_produced == total_out)
                    wfsm <= W_IDLE;
            end

            default: wfsm <= W_IDLE;
            endcase
        end
    end

    // ── Squeeze side: unpack kt_dout (32b) into 4 sequential bytes ───────────
    reg [31:0] sq_word;
    reg [1:0]  sq_byte_idx;
    reg        sq_word_valid;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            sq_word_valid <= 0;
            sq_byte_idx   <= 0;
        end else if (wfsm == W_SQUEEZE) begin
            if (kt_dout_valid && kt_dout_ready) begin
                sq_word       <= kt_dout;
                sq_word_valid <= 1;
                sq_byte_idx   <= 0;
            end else if (sq_word_valid && squeeze_ready && (out_produced < total_out)) begin
                if (sq_byte_idx == 2'd3) sq_word_valid <= 0;
                else                     sq_byte_idx   <= sq_byte_idx + 2'd1;
            end

            if (sq_word_valid && squeeze_ready && (out_produced < total_out))
                out_produced <= out_produced + 1'b1;
        end else begin
            sq_word_valid <= 0;
            sq_byte_idx   <= 0;
        end
    end

    assign kt_dout_ready  = (wfsm == W_SQUEEZE) && !sq_word_valid && (out_produced < total_out);
    assign squeeze_valid  = (wfsm == W_SQUEEZE) && sq_word_valid && (out_produced < total_out);
    assign squeeze_data   = sq_word[sq_byte_idx*8 +: 8];

endmodule
