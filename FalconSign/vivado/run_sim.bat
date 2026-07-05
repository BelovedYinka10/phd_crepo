@echo off
:: run_sim.bat — FalconSign SamplerZ simulation with Vivado XSim on Windows
::
:: Usage (from Vivado 2023.2 Tcl Shell or after sourcing settings):
::   cd FalconSign\vivado
::   run_sim.bat
::
:: Or open "Vivado 2023.2 Tcl Shell" from the Start menu, navigate here, then run.
::
:: What it does:
::   1. Generates memory.hex (needs Python 3 in PATH)
::   2. Compiles all SystemVerilog with xvlog
::   3. Elaborates with xelab
::   4. Runs simulation with xsim (prints pass/fail to console)

setlocal

set ROOT=..
set INC=%ROOT%\include
set SAM=%ROOT%\sampler
set SIM=%ROOT%\sampler_sim
set TOP=samplerz_tb
set SNAP=samplerz_tb_snap

:: ── Step 0: generate memory.hex ─────────────────────────────────────────────
echo === Generating memory.hex ===
python %SIM%\gen_mem.py memory.hex
if errorlevel 1 (
    echo ERROR: Failed to generate memory.hex
    echo Make sure Python 3 is installed and in PATH.
    exit /b 1
)

:: ── Step 1: compile all SystemVerilog ───────────────────────────────────────
:: All files in ONE xvlog call so include guards work correctly.
:: Package files (falconsoar_pkg, sample_pkg) must come first.
echo === Compiling SystemVerilog sources ===
xvlog -sv -i %INC% ^
    %INC%\falconsoar_pkg.sv ^
    %INC%\sample_pkg.sv ^
    %INC%\time_counter.sv ^
    %INC%\fpr_mul_opt_2pipe.sv ^
    %INC%\fpr_mul_const_opt_2pipe.sv ^
    %INC%\floating_point_add_2pip.sv ^
    %INC%\floating_point_sub_2pip.sv ^
    %INC%\fp_sub_s.sv ^
    %INC%\fp_mult_s.sv ^
    %INC%\fp_flt2i_int32_s.sv ^
    %INC%\fp_i2flt_int32_s.sv ^
    %INC%\fp_i2flt_int32.sv ^
    %INC%\fp_flt2i_int64_s.sv ^
    %INC%\fp_i2flt_int64_s.sv ^
    %INC%\counter_ce.sv ^
    %INC%\rbsh.sv ^
    %SAM%\chacha20.sv ^
    %SAM%\berexp.sv ^
    %SAM%\fpr_cal.sv ^
    %SAM%\pre_samp.sv ^
    %SAM%\samp_loop.sv ^
    %SAM%\refill_control.sv ^
    %SAM%\samplerz.sv ^
    %SAM%\samplerz_top.sv ^
    %SIM%\samplerz_tb.sv

if errorlevel 1 (
    echo ERROR: xvlog compilation failed. Check output above.
    exit /b 1
)

:: ── Step 2: elaborate ────────────────────────────────────────────────────────
echo === Elaborating: top=%TOP% ===
xelab -debug typical %TOP% -s %SNAP%
if errorlevel 1 (
    echo ERROR: xelab elaboration failed.
    exit /b 1
)

:: ── Step 3: simulate ─────────────────────────────────────────────────────────
echo === Running simulation ===
xsim %SNAP% --runall
:: xsim returns non-zero on $finish — that is normal, not an error
echo === Simulation complete ===

endlocal
