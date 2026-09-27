# Radar Correlation & DSP Suite

This repository implements a high-precision, hardware-efficient Radar processing chain in standard VHDL. It features a bit-true workflow validating DSP48E1 hardware implementations against Octave/MATLAB golden reference models.

The primary engine is a **Dual-Channel Quadrature (I/Q) Range Detector / Matched Filter** utilizing two parallel 5x time-division folded 50-tap symmetric FIR architectures. Running at 250 MHz on budget Zynq-7000 silicon, the complete dual-channel processing chain requires only **10 DSP48E1 slices** while maintaining full timing closure.

---

## Repository Structure

```text
Correlation/
├── constraints/                  # Timing and physical placement constraints
│   └── range_detector.xdc        # 250 MHz clock definition and CDC path constraints
├── docs/                         # Hardware verification snapshots and documentation
│   ├── Firt_matched_filter.png   # Single-channel matched filter waveform snapshot
│   ├── IQ_match_filter.png       # Dual-channel (I/Q) matched filter waveform snapshot
│   ├── tx_noisy_rx_pulses.png    # Octave simulation: Transmitted and noisy RX pulses
│   └── detected_range.png        # Octave simulation: Coherent integration range detection
├── ipcores/                      # AMD/Xilinx IP core definitions
│   ├── adc_fifo/                 # Dual-clock asynchronous CDC FIFO
│   ├── clk_dsp/                  # Clocking Wizard generating 250 MHz DSP clock
│   └── IQ2vector/                # CORDIC phase extraction IP core
├── m_files/                      # Octave/MATLAB scripts for bit-true golden models
│   ├── Correlation.m             # System-level correlation & matched filter simulation model
│   ├── golden_ref_vector.dat     # Exported golden reference test data
│   ├── my_comp_round.m           # Complex rounding helper model
│   ├── myconv.m                  # Bit-true fixed-point convolution algorithm
│   └── myround.m                 # Golden reference for convergent (round-to-even) rounding
├── rtl_sim/                      # Testbenches and simulation assets
│   ├── IQ2Vect_tb.vhd            # Vectoring module testbench
│   ├── range_detector.wcfg       # Vivado XSIM waveform configuration
│   ├── rtl_golden_ref_vector.vhd # Golden reference test vectors
│   ├── tb_range_detector.vhd     # Top-level range detector testbench
│   └── test_fix_round.vhd        # Fixed-point rounding testbench
├── rtl_src/                      # Synthesizable RTL source files
│   ├── conv_rounding.vhd         # Convergent rounding (round-to-even) logic
│   ├── DSP_wrapper.vhd           # Parametric DSP48E1 macro wrapper
│   ├── fir_impl.vhd              # 50-tap symmetric folded FIR filter engine
│   ├── radar_pkg.vhd             # Design constants, types, and FIR coefficients
│   ├── range_detector.vhd        # Dual-channel CDC wrapper and top-level processing
│   └── vectoring.vhd             # CORDIC magnitude and phase processing interface
├── scripts/                      # Project build automation
│   ├── range_detector.tcl        # Vivado Project Mode generation script
│   └── runme.bat                 # One-click Windows batch installer
├── .gitignore                    # Standard Git exclusions for Vivado build artifacts
├── LICENSE                       # Project licensing terms
└── README.md                     # Repository documentation
```

---

## MATLAB / Octave System Modeling (`m_files`)

The `m_files/` directory contains system-level scripts used to design, quantize, and verify the VHDL RTL implementation:

1. **System Correlation Model (`Correlation.m`):** Top-level simulation modeling signal pulse compression, matched filtering, and theoretical correlation limits.
2. **Convolution & Fixed-Point DSP (`myconv.m`):** Custom convolution routine modeling fixed-point arithmetic before RTL migration.
3. **Bit-True Rounding Models (`myround.m` & `my_comp_round.m`):** Implements convergent rounding (round-to-even) and complex rounding models to eliminate DC bias during bit-width reduction.
4. **Golden Vector Dataset (`golden_ref_vector.dat`):** Exported test vector dataset used to populate `rtl_sim/rtl_golden_ref_vector.vhd` for self-checking VHDL simulation.

### Simulation & Golden Vector Workflow
The correlation simulation is controlled via `Correlation(act_prnt, new_test_vector)` in Octave/MATLAB:

* **Debugging & Printing (`act_prnt`):** Set to `1` to enable verbose console debugging and dump VHDL-ready ROM arrays (`chrp_ampl_c`, `chrp_phs_c`, `inph_c`, `quadr_c`).
* **Vector Generation (`new_test_vector`):**
  * `Correlation(0, 1)`: Generates a new stochastic noisy RX signal package and writes updated reference files (`golden_ref_vector.m` for Octave and `rtl_golden_ref_vector.vhd` for the Vivado RTL testbench).
  * `Correlation(0, 0)`: Loads the static baseline vector set. Running this mode guarantees that Octave and RTL simulations operate on an identical, unified test dataset for strict bit-true correlation verification.

> **Note:** A pre-generated golden vector file (`rtl_sim/rtl_golden_ref_vector.vhd`) is provided in the repository as the default baseline reference.

### Simulation & Theoretical Limits
The following snapshots demonstrate the Octave simulation results, establishing the theoretical baseline for the RTL implementation:

![Transmitted and Noisy RX Pulses](docs/tx_noisy_rx_pulses.png)
*Transmitted chirp and deeply embedded noisy RX signal (SNR = -30.25dB).*

![Coherent Integration and Range Detection](docs/detected_range.png)
*Matched filter output showcasing successful pulse compression and coherent integration gain.*

---

## Hardware Architecture & Verification

### 1. Dual-Channel Quadrature (I/Q) Matched Filter
* **Parallel Processing Engines:** Dual 50-tap symmetric FIR filters process In-Phase (I) and Quadrature (Q) ADC channels in parallel.
* **5x Time-Division Folding:** Each channel folds 50 taps down to 5 cascaded DSP48E1 slices operating at 250 MHz (5 clock cycles per input sample).
* **Resource Optimization:** Consumes only 10 DSP48E1 slices (~4.5% of Zynq-7000 DSP resources) for complete complex I/Q filtering while achieving zero timing violations.
* **Clock Domain Crossing (CDC):** Dual ADC channels cross safely into the 250 MHz DSP clock domain via an asynchronous FIFO (`adc_fifo`).

### 2. Waveform Verification

#### Single-Channel Baseline Verification
Functional baseline validation for the single-channel folded FIR engine:

![Single-Channel Matched Filter Waveform Output](docs/Firt_matched_filter.png)

#### Dual-Channel Quadrature (I/Q) Verification
Parallel channel output validation for In-Phase and Quadrature paths (`fir_amp_s` and `fir_pha_s`) against the Octave golden reference vector package (`rtl_golden_ref_vector.vhd`):

![Quadrature Matched Filter Waveform Output](docs/IQ_match_filter.png)

---

## Build Instructions

To generate the complete Vivado project environment using the automated project mode script:

1. Open a command prompt inside `scripts/`.
2. Run the batch launcher:

```cmd
runme.bat
```

The script automatically:
* Sources the Xilinx Vivado environment via `settings64.bat`.
* Imports IP cores (`adc_fifo.xci`, `clk_dsp.xci`).
* Adds targeted RTL sources from `rtl_src/` and simulation assets from `rtl_sim/`.
* Sets `range_detector` as the top-level entity and `tb_range_detector` as the simulation top.
* Applies constraints from `constraints/range_detector.xdc`.