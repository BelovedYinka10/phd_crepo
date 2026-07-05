// counter_ce.sv
// Parameterized up-counter with clock enable and configurable terminal count.
//
// Counts from 0 to num-1 (num cycles total), then pulses done and stops.
// When num == 0 the counter wraps naturally, producing exactly 2^BW counts.
// This is intentional: callers pass num = BW'(2^BW) which truncates to 0.
//
// Port summary
//   start  : synchronous start; resets cnt to 0 and begins counting
//   num    : terminal count (exclusive); 0 means 2^BW
//   run    : clock enable; counting pauses while low
//   cnt    : current counter value
//   cnt_en : high for every cycle that cnt holds a valid output
//   last   : combinational; high on the last valid count (cnt == num-1)
//   done   : registered one-cycle pulse, one cycle after last

`timescale 1ns/1ps

module counter_ce #(
    parameter int BW = 8
) (
    input  wire          clk,
    input  wire          rst_n,
    input  wire          start,
    input  wire [BW-1:0] num,      // 0 = count 2^BW cycles
    input  wire          run,      // clock enable
    output reg  [BW-1:0] cnt,
    output reg           cnt_en,   // high while counting
    output wire          last,     // high on final count cycle
    output reg           done      // one-cycle pulse after final count
);

    // num=0 → end_val = all-ones = 2^BW-1 (natural wrap)
    wire [BW-1:0] end_val = num - 1'b1;

    assign last = cnt_en & run & (cnt == end_val);

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            cnt    <= '0;
            cnt_en <= 1'b0;
            done   <= 1'b0;
        end else begin
            done <= last;   // registered: fires one cycle after the last valid cnt
            if (start) begin
                cnt    <= '0;
                cnt_en <= 1'b1;
            end else if (last) begin
                cnt    <= '0;
                cnt_en <= 1'b0;
            end else if (cnt_en && run) begin
                cnt    <= cnt + 1'b1;
            end
        end
    end

endmodule