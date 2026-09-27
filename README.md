# Radar Correlation & DSP Suite

This repository implements a high-precision, hardware-efficient Radar processing chain in standard VHDL. It features a bit-true workflow validating DSP48E1 hardware implementations against Octave/MATLAB golden reference models.

The primary engine is a **Range Detector / Matched Filter** using a 5x time-division folded 50-tap symmetric FIR architecture. Running at 250 MHz on budget Zynq-7000 silicon, the engine requires only **5 DSP48E1 slices** while achieving full timing closure.

---

## Repository Structure

```text
Correlation/
├── constraints/                  # Timing and physical placement constraints
│   └── range_detector.xdc        # 250 MHz clock definition and CDC path constraints
├── docs/                         # Hardware verification snapshots and documentation
│   └── Firt_matched_filter.png   # Simulation waveform snapshot of matched filter output
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
│   ├── fir_impl.vhd              # 50-tap symmetric folded FIR filter (5 DSP slices)
│   ├── radar_pkg.vhd             # Design constants, types, and FIR coefficients
│   ├── range_detector.vhd        # Top-level CDC wrapper and threshold detector
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

The `m_files/` directory contains system-level scripts[cite: 14] used to design, quantize, and verify the VHDL RTL implementation:

1. **System Correlation Model (`Correlation.m`):** Top-level simulation modeling signal pulse compression, matched filtering, and theoretical correlation limits[cite: 14].
2. **Convolution & Fixed-Point DSP (`myconv.m`):** Custom convolution routine modeling fixed-point arithmetic before RTL migration[cite: 14].
3. **Bit-True Rounding Models (`myround.m` & `my_comp_round.m`):** Implements convergent rounding (round-to-even) and complex rounding models[cite: 14] to eliminate DC bias during bit-width reduction.
4. **Golden Vector Dataset (`golden_ref_vector.dat`):** Exported test vector dataset[cite: 14] used to populate `rtl_sim/rtl_golden_ref_vector.vhd` for self-checking VHDL simulation.

---

## Hardware Architecture & Verification

### 1. 50-Tap Folded Matched Filter
* **5x Time-Division Folding:** Reduces a 50-tap symmetric FIR filter down to 5 cascaded DSP48E1 slices operating at 250 MHz (5 clock cycles per input sample).
* **Resource Optimization:** Consumes only 5 DSP48E1 slices (~2.27% of Zynq-7000 DSP resources), 414 LUTs, and 924 FFs while achieving zero timing violations.
* **Clock Domain Crossing (CDC):** Input ADC samples cross safely into the 250 MHz DSP clock domain via an asynchronous FIFO (`adc_fifo`).

### 2. Waveform Verification
Functional accuracy is validated against the Octave golden reference vector package (`rtl_golden_ref_vector.vhd`) in simulation.

![Matched Filter Waveform Output](docs/Firt_matched_filter.png)

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