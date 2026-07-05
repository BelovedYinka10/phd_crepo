// shake256.sv — SHAKE256 extendable-output function (NIST FIPS 202)
//
// Rate = 136 bytes (1088 bits), capacity = 512 bits.
// Domain separation suffix: 0x1F  (SHAKE256, FIPS 202 §6.3).
//
// Instantiates keccak_f1600 for the permutation (24 cycles/call).
//
// Protocol:
//   1. Pulse init to begin a new hash (also effective after rst_n).
//   2. Feed message bytes on absorb_data while asserting absorb_valid.
//      Set absorb_last = 1 on the final byte.
//      The module consumes a byte only when ready = 1.
//   3. Wait for busy = 0.
//   4. Pulse squeeze_ready once per desired output byte.
//      Read squeeze_data when squeeze_valid = 1 (one cycle later).
//      The XOF re-permutes automatically every 136 output bytes.
//   5. To hash a new message, pulse init before step 2.
//
// Timing: perm_start is a 1-cycle pulse; perm_done appears 25 cycles after
// perm_start (24 rounds, each 1 cycle, state registered on entry).

`timescale 1ns/1ps

module shake256 (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        init,          // 1-cycle pulse: clear state, begin new hash
    // absorb
    input  wire        absorb_valid,
    input  wire [7:0]  absorb_data,
    input  wire        absorb_last,   // high on final message byte
    output wire        ready,         // 1 = accepting message bytes
    // squeeze
    input  wire        squeeze_ready, // pulse: request next output byte
    output reg  [7:0]  squeeze_data,
    output reg         squeeze_valid, // 1 cycle after squeeze_ready: data valid
    // status
    output wire        busy           // 1 while Keccak permutation running
);

    localparam RATE = 136;  // 1088 bits / 8

    // ── Keccak-f[1600] ────────────────────────────────────────────────────────
    reg  [1599:0] state;
    wire [1599:0] perm_out;
    reg           perm_start;
    wire          perm_done;

    keccak_f1600 u_keccak (
        .clk       (clk),
        .rst_n     (rst_n),
        .start     (perm_start),
        .state_in  (state),
        .state_out (perm_out),
        .done      (perm_done)
    );

    // ── FSM ───────────────────────────────────────────────────────────────────
    localparam [1:0] S_ABSORB  = 2'd0;
    localparam [1:0] S_PERMUTE = 2'd1;
    localparam [1:0] S_SQUEEZE = 2'd2;

    reg [1:0]  fsm;
    reg [1:0]  after_perm;  // FSM state to enter once perm_done fires
    reg [7:0]  byte_pos;    // absorb: current byte in rate block  (0 .. RATE-1)
    reg [7:0]  sq_pos;      // squeeze: current byte in rate block (0 .. RATE-1)

    assign ready = (fsm == S_ABSORB);
    assign busy  = (fsm == S_PERMUTE);

    // ── State byte helpers ────────────────────────────────────────────────────
    // FIPS 202 §3.1: byte b of state is in lane b/8 at bit offset (b%8)*8.
    // Lane i occupies state[i*64 +: 64].

    function automatic [1599:0] xor_byte;
        input [1599:0] st;
        input [7:0]    pos;
        input [7:0]    val;
        integer lane, boff;
        begin
            lane  = pos / 8;
            boff  = (pos % 8) * 8;
            xor_byte = st;
            xor_byte[lane*64 + boff +: 8] = st[lane*64 + boff +: 8] ^ val;
        end
    endfunction

    function automatic [7:0] read_byte;
        input [1599:0] st;
        input [7:0]    pos;
        integer lane, boff;
        begin
            lane      = pos / 8;
            boff      = (pos % 8) * 8;
            read_byte = st[lane*64 + boff +: 8];
        end
    endfunction

    // ── Sequential logic ──────────────────────────────────────────────────────
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state         <= '0;
            fsm           <= S_ABSORB;
            after_perm    <= S_ABSORB;
            byte_pos      <= '0;
            sq_pos        <= '0;
            perm_start    <= 0;
            squeeze_data  <= '0;
            squeeze_valid <= 0;
        end else begin
            perm_start    <= 0;      // default: deasserted every cycle
            squeeze_valid <= 0;

            if (init) begin
                // Reset for a new message. If a permutation is in flight it will
                // complete but its result is discarded (perm_done ignored in S_ABSORB).
                state    <= '0;
                byte_pos <= '0;
                sq_pos   <= '0;
                fsm      <= S_ABSORB;
            end else begin
                case (fsm)

                // ── Absorb ────────────────────────────────────────────────────
                S_ABSORB: begin
                    if (absorb_valid) begin
                        if (absorb_last) begin
                            // Last byte: absorb it, append SHAKE256 suffix 0x1F at
                            // byte_pos, append 0x80 at RATE-1 (multi-rate padding).
                            // Both XORs are combined into a single register write.
                            state      <= xor_byte(
                                              xor_byte(state, byte_pos,
                                                       absorb_data ^ 8'h1F),
                                              RATE - 1, 8'h80);
                            after_perm <= S_SQUEEZE;
                            perm_start <= 1;
                            byte_pos   <= '0;
                            fsm        <= S_PERMUTE;
                        end else if (byte_pos == RATE - 1) begin
                            // Rate block full mid-message: permute, then resume absorb.
                            state      <= xor_byte(state, byte_pos, absorb_data);
                            after_perm <= S_ABSORB;
                            perm_start <= 1;
                            byte_pos   <= '0;
                            fsm        <= S_PERMUTE;
                        end else begin
                            state    <= xor_byte(state, byte_pos, absorb_data);
                            byte_pos <= byte_pos + 1;
                        end
                    end
                end

                // ── Wait for Keccak-f[1600] ───────────────────────────────────
                S_PERMUTE: begin
                    if (perm_done) begin
                        state  <= perm_out;
                        sq_pos <= '0;        // reset on every permutation completion
                        fsm    <= after_perm;
                    end
                end

                // ── Squeeze ───────────────────────────────────────────────────
                S_SQUEEZE: begin
                    if (squeeze_ready) begin
                        squeeze_data  <= read_byte(state, sq_pos);
                        squeeze_valid <= 1;
                        if (sq_pos == RATE - 1) begin
                            // Output block exhausted: permute for next 136 bytes.
                            after_perm <= S_SQUEEZE;
                            perm_start <= 1;
                            fsm        <= S_PERMUTE;
                            sq_pos     <= '0;
                        end else begin
                            sq_pos <= sq_pos + 1;
                        end
                    end
                end

                default: fsm <= S_ABSORB;
                endcase
            end
        end
    end

endmodule
