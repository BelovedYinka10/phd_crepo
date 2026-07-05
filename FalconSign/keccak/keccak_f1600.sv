// keccak_f1600.sv — Keccak-f[1600] permutation (NIST FIPS 202)
//
// Applies 24 rounds, one round per clock cycle (24 cycles total).
// State: 1600 bits = 25 × 64-bit lanes, lane[x+5*y] at bits [64*(x+5*y)+63:64*(x+5*y)].
// Byte ordering: lane[i] bit 0 = byte 8*i bit 0 (little-endian lanes, per FIPS 202 §3.1).
//
// Interface:
//   start  — 1-cycle pulse: load state_in, begin 24 rounds
//   done   — 1-cycle pulse on completion (24 cycles after start)
//   state_out — valid when done is high, stable until next start

`timescale 1ns/1ps

module keccak_f1600 (
    input  wire          clk,
    input  wire          rst_n,
    input  wire          start,
    input  wire [1599:0] state_in,
    output reg  [1599:0] state_out,
    output reg           done
);

    // ── Round constants (FIPS 202 Table 1) ───────────────────────────────────
    localparam [63:0] RC [0:23] = '{
        64'h0000000000000001, 64'h0000000000008082,
        64'h800000000000808A, 64'h8000000080008000,
        64'h000000000000808B, 64'h0000000080000001,
        64'h8000000080008081, 64'h8000000000008009,
        64'h000000000000008A, 64'h0000000000000088,
        64'h0000000080008009, 64'h000000008000000A,
        64'h000000008000808B, 64'h800000000000008B,
        64'h8000000000008089, 64'h8000000000008003,
        64'h8000000000008002, 64'h8000000000000080,
        64'h000000000000800A, 64'h800000008000000A,
        64'h8000000080008081, 64'h8000000000008080,
        64'h0000000080000001, 64'h8000000080008008
    };

    // ── State ─────────────────────────────────────────────────────────────────
    reg [63:0] s [0:24];
    reg [4:0]  round;
    reg        running;

    // ── Left-rotate helper ────────────────────────────────────────────────────
    function automatic [63:0] rotl;
        input [63:0] v;
        input [5:0]  n;
        rotl = (n == 0) ? v : {v[63-n:0], v[63:64-n]};
    endfunction

    // ── One-round combinational logic ─────────────────────────────────────────
    function automatic void keccak_round(
        input  [63:0] A  [0:24],
        input  [63:0] rc,
        output [63:0] A2 [0:24]
    );
        reg [63:0] C[0:4], D[0:4], B[0:24];
        // θ
        for (int x = 0; x < 5; x++)
            C[x] = A[x] ^ A[x+5] ^ A[x+10] ^ A[x+15] ^ A[x+20];
        for (int x = 0; x < 5; x++)
            D[x] = C[(x+4)%5] ^ rotl(C[(x+1)%5], 1);
        for (int i = 0; i < 25; i++)
            B[i] = A[i] ^ D[i%5];
        // ρ + π  (pre-computed mapping: dest ← rotl(B[src], rot))
        A2[ 0] = rotl(B[ 0],  0); A2[ 1] = rotl(B[ 6], 44);
        A2[ 2] = rotl(B[12], 43); A2[ 3] = rotl(B[18], 21);
        A2[ 4] = rotl(B[24], 14); A2[ 5] = rotl(B[ 3], 28);
        A2[ 6] = rotl(B[ 9], 20); A2[ 7] = rotl(B[10],  3);
        A2[ 8] = rotl(B[16], 45); A2[ 9] = rotl(B[22], 61);
        A2[10] = rotl(B[ 1],  1); A2[11] = rotl(B[ 7],  6);
        A2[12] = rotl(B[13], 25); A2[13] = rotl(B[19],  8);
        A2[14] = rotl(B[20], 18); A2[15] = rotl(B[ 4], 27);
        A2[16] = rotl(B[ 5], 36); A2[17] = rotl(B[11], 10);
        A2[18] = rotl(B[17], 15); A2[19] = rotl(B[23], 56);
        A2[20] = rotl(B[ 2], 62); A2[21] = rotl(B[ 8], 55);
        A2[22] = rotl(B[14], 39); A2[23] = rotl(B[15], 41);
        A2[24] = rotl(B[21],  2);
        // χ — reads must come from pre-χ state; copy before modifying in place
        begin
            reg [63:0] T [0:24];
            for (int i = 0; i < 25; i++) T[i] = A2[i];
            for (int x = 0; x < 5; x++)
                for (int y = 0; y < 5; y++)
                    A2[x+5*y] = T[x+5*y] ^ ((~T[(x+1)%5+5*y]) & T[(x+2)%5+5*y]);
        end
        // ι
        A2[0] = A2[0] ^ rc;
    endfunction

    // ── FSM ───────────────────────────────────────────────────────────────────
    reg [63:0] s_next [0:24];

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            running <= 0; round <= 0; done <= 0;
        end else begin
            done <= 0;
            if (start) begin
                for (int i = 0; i < 25; i++) s[i] <= state_in[64*i+:64];
                round   <= 0;
                running <= 1;
            end else if (running) begin
                keccak_round(s, RC[round], s_next);
                for (int i = 0; i < 25; i++) s[i] <= s_next[i];
                if (round == 5'd23) begin
                    running <= 0;
                    done    <= 1;
                end else begin
                    round <= round + 1;
                end
            end
        end
    end

    always_comb begin
        for (int i = 0; i < 25; i++) state_out[64*i+:64] = s[i];
    end

endmodule
