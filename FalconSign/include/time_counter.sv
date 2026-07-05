// time_counter.sv
// Simple cycle counter: counts clock cycles from start to done.
// Used in samplerz_top for performance measurement (not functionally critical).

`timescale 1ns/1ps

module time_counter (
    input  wire clk,
    input  wire rst_n,
    input  wire start,
    input  wire done
);
    reg [31:0] cnt;
    reg        running;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            cnt     <= '0;
            running <= 1'b0;
        end else if (start) begin
            cnt     <= '0;
            running <= 1'b1;
        end else if (done) begin
            running <= 1'b0;
        end else if (running) begin
            cnt <= cnt + 1'b1;
        end
    end

endmodule
