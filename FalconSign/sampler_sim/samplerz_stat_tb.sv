`timescale 1ns/1ps
// samplerz_stat_tb.sv — Statistical correctness testbench.
//
// Runs N_TASKS pairs of samples (each task produces 2 samples),
// prints each result as:   SAMPLE: <128-bit hex word>
// where bits [63:0]   = first  sample (IEEE 754 double of integer z0)
//       bits [127:64] = second sample (IEEE 754 double of integer z1)
//
// stat_test.py reads this output and runs a chi-squared test against
// the Python samplerz reference implementation.

`include "falconsoar_pkg.sv"
`include "sample_pkg.sv"

module samplerz_stat_tb
    import falconsoar_pkg::*;
    import sample_pkg::*;
;
    parameter int N_TASKS = 1500;  // → 3000 samples total

    logic clk   = 0;
    logic rst_n = 0;
    always #5 clk = ~clk;

    exec_operator_if task_itf ();
    mem_inst_if      mem_rd   ();
    mem_inst_if      mem_wr   ();

    samplerz_top dut (
        .clk    (clk),
        .rst_n  (rst_n),
        .task_itf(task_itf),
        .mem_rd  (mem_rd),
        .mem_wr  (mem_wr)
    );

    localparam int MEM_WORDS = BANK_NUM * BANK_DEPTH;

    logic [BANK_WIDTH-1:0] mem_array [MEM_WORDS];

    // 4-cycle read latency matching SAMPLERZ_READ_DELAY=4.
    // Address is sampled unconditionally every cycle (no enable gating) to avoid
    // X-propagation from pulse_extender flip-flops that have no reset.
    // BRAM-like: output holds last read value.
    logic [BANK_WIDTH-1:0] rd_pipe [4];

    always_ff @(posedge clk) begin
        rd_pipe[0] <= mem_array[mem_rd.addr];
        rd_pipe[1] <= rd_pipe[0];
        rd_pipe[2] <= rd_pipe[1];
        rd_pipe[3] <= rd_pipe[2];
    end
    assign mem_rd.data = rd_pipe[3];

    always_ff @(posedge clk)
        if (mem_wr.en) mem_array[mem_wr.addr] <= mem_wr.data;

    localparam [MEM_ADDR_BITS-1:0] SRC0 = 13'd0;
    localparam [MEM_ADDR_BITS-1:0] SRC1 = 13'd2;
    localparam [MEM_ADDR_BITS-1:0] DST  = 13'd64;

    localparam [TASK_REDUCE_BW-1:0] TASK_INIT =
        {SRC0, SRC1, DST, 13'd0, 1'b1, 4'd1, 11'd0};
    localparam [TASK_REDUCE_BW-1:0] TASK_SAMP =
        {SRC0, SRC1, DST, 13'd0, 1'b0, 4'd1, 11'd0};

    task tick; @(posedge clk); #1; endtask

    integer timeout, i;

    initial begin
        task_itf.start      = 1'b0;
        task_itf.input_task = '0;
        $readmemh("memory.hex", mem_array);

        repeat(10) tick;
        rst_n = 1'b1;
        repeat(5) tick;

        // ── INIT: start ChaCha20 PRNG ────────────────────────────────────────
        task_itf.start      = 1'b1;
        task_itf.input_task = TASK_INIT;
        tick;
        task_itf.start = 1'b0;
        timeout = 0;
        while (!task_itf.op_done && timeout < 200_000) begin tick; timeout++; end
        if (!task_itf.op_done) begin $display("ERROR: INIT timeout"); $finish(1); end
        tick; tick;

        // ── Run N_TASKS, each produces 2 samples ────────────────────────────
        for (i = 0; i < N_TASKS; i++) begin
            task_itf.start      = 1'b1;
            task_itf.input_task = TASK_SAMP;
            tick;
            task_itf.start = 1'b0;

            timeout = 0;
            while (!task_itf.op_done && timeout < 200_000) begin tick; timeout++; end
            if (!task_itf.op_done) begin
                $display("ERROR: task %0d timeout", i);
                $finish(1);
            end
            tick;  // let memory write settle

            // Print: bits[63:0]=first_sample, bits[127:64]=second_sample
            $display("SAMPLE: %032X", mem_array[DST][127:0]);
        end

        $display("DONE: %0d tasks (%0d samples)", N_TASKS, N_TASKS*2);
        $finish(0);
    end

endmodule
