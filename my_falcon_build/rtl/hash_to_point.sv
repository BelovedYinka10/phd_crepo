// hash_to_point.sv — Falcon hash_to_point (FIPS 202 / Falcon spec §3.7.1)
//
// Converts (salt || message) via SHAKE256 into n coefficients in Z_q = Z/12289.
// Rejection sampling: read 16-bit BE words from SHAKE256 XOF;
//   accept if word < k*Q  (k=5, Q=12289, threshold=61445);
//   output = accepted_word mod Q.
//
// Interface:
//   n_log2 [1:0]  — 0→n=512 (Falcon-512), 1→n=1024 (Falcon-1024)
//   msg_valid / msg_data / msg_last — byte stream (salt then message)
//   coeff_valid / coeff_data [13:0] — output coefficient stream  (one per cycle)
//   done           — 1-cycle pulse when all n coefficients delivered
//   busy           — high while hashing

`timescale 1ns/1ps

module hash_to_point (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        start,        // 1-cycle pulse: begin new hash
    input  wire [1:0]  n_log2,       // 0=n512, 1=n1024
    // message input (absorb phase)
    input  wire        msg_valid,
    input  wire [7:0]  msg_data,
    input  wire        msg_last,     // pulse on final byte of (salt||message)
    // coefficient output (squeeze phase)
    output reg         coeff_valid,
    output reg  [13:0] coeff_data,   // value in [0, Q-1]
    output reg         done,
    output wire        busy,
    output wire        msg_ready     // =1 when accepting message bytes
);
    localparam [13:0] Q         = 14'd12289;
    localparam [16:0] THRESHOLD = 17'd61445;  // k*Q = 5*12289

    // ── SHAKE256 instance ─────────────────────────────────────────────────────
    wire        sh_absorb_valid, sh_absorb_last, sh_sq_ready;
    wire [7:0]  sh_absorb_data;
    wire [7:0]  sh_sq_data;
    wire        sh_sq_valid, sh_busy, sh_ready;

    shake256 u_shake (
        .clk          (clk),
        .rst_n        (rst_n),
        .absorb_valid (sh_absorb_valid),
        .absorb_data  (sh_absorb_data),
        .absorb_last  (sh_absorb_last),
        .squeeze_ready(sh_sq_ready),
        .squeeze_data (sh_sq_data),
        .squeeze_valid(sh_sq_valid),
        .busy         (sh_busy),
        .ready        (sh_ready)
    );

    assign sh_absorb_valid = msg_valid && sh_ready;
    assign sh_absorb_data  = msg_data;
    assign sh_absorb_last  = msg_last;
    assign msg_ready       = sh_ready;

    // ── Squeeze + rejection sampling FSM ─────────────────────────────────────
    localparam S_IDLE    = 2'd0;
    localparam S_ABSORB  = 2'd1;
    localparam S_SQUEEZE = 2'd2;
    localparam S_DONE    = 2'd3;

    reg [1:0]  fsm;
    reg [10:0] coeff_cnt;   // how many coefficients accepted so far
    reg [10:0] n_target;    // n (512 or 1024)
    reg        hi_byte_valid;
    reg [7:0]  hi_byte;
    reg        sq_req;      // request next byte from SHAKE256

    assign busy = (fsm != S_IDLE);
    assign sh_sq_ready = sq_req;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            fsm           <= S_IDLE;
            coeff_cnt     <= '0;
            hi_byte_valid <= 0;
            coeff_valid   <= 0;
            done          <= 0;
            sq_req        <= 0;
        end else begin
            coeff_valid <= 0;
            done        <= 0;
            sq_req      <= 0;

            case (fsm)
            S_IDLE: begin
                if (start) begin
                    n_target      <= (n_log2 == 0) ? 11'd512 : 11'd1024;
                    coeff_cnt     <= '0;
                    hi_byte_valid <= 0;
                    fsm           <= S_ABSORB;
                end
            end

            // Wait for message to be fully absorbed
            S_ABSORB: begin
                if (!sh_busy && !sh_ready) begin
                    // SHAKE256 is in squeeze state
                    fsm    <= S_SQUEEZE;
                    sq_req <= 1;   // request first byte
                end
            end

            S_SQUEEZE: begin
                if (sh_sq_valid) begin
                    if (!hi_byte_valid) begin
                        // First byte of a 16-bit word
                        hi_byte       <= sh_sq_data;
                        hi_byte_valid <= 1;
                        sq_req        <= 1;  // request low byte
                    end else begin
                        // Second byte — complete the 16-bit word
                        automatic reg [16:0] word;
                        word = {1'b0, hi_byte, sh_sq_data}; // big-endian
                        hi_byte_valid <= 0;
                        if (word < THRESHOLD) begin
                            // Accept: output coefficient mod Q
                            coeff_data  <= word % Q;  // word < 5*Q < 2^14 after mod
                            coeff_valid <= 1;
                            coeff_cnt   <= coeff_cnt + 1;
                            if (coeff_cnt + 1 == n_target) begin
                                fsm  <= S_DONE;
                                done <= 1;
                            end else begin
                                sq_req <= 1;  // request next byte
                            end
                        end else begin
                            sq_req <= 1;  // reject: try again
                        end
                    end
                end
            end

            S_DONE: fsm <= S_IDLE;
            endcase
        end
    end

endmodule
