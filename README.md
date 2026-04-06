# Radar Correlation & DSP Suite

This project implements a high-precision Radar processing chain, bridging the gap between theoretical DSP and hardware-efficient RTL. It features a **Bit-True workflow** where Octave models are used to validate VHDL implementations of Matched Filters and CORDIC pre-processing components.

## 📂 Repository Structure

* **`m_files/`**: Octave/MATLAB scripts for system-level modeling, including the `myround` golden reference for bit-true hardware verification.
* **`RTL_Src/`**: Synthesizable VHDL source code. Contains the core logic, including the Fixed-Point Rounding and Symmetric Saturation modules.
* **`RTL_Sim/`**: VHDL testbenches, simulation scripts, and `xsim` waveform configuration files for functional verification.
* **`Constraint/`**: FPGA constraint files (.xdc) for defining clock timing and physical pin assignments.
* **`.gitignore`**: Standard Git configuration to exclude simulation artifacts and Vivado temporary files.
* **`LICENSE`**: Project licensing terms.

## 🎯 Project Goals & Features

The objective is to implement a robust Radar processor capable of high dynamic range and accurate phase detection. Key technical highlights include:

1.  **Pulse Compression**: Matched filtering of LFM (Chirp) waveforms to maximize SNR.
2.  **Bit-True Arithmetic**: Implementation of **Convergent Rounding (Round-to-Even)** to eliminate DC bias during bit-width reduction (31-bit to 16-bit).
3.  **Dynamic Range Management**: Hardware-level **Symmetric Saturation** logic to prevent Two's Complement wrap-around and catastrophic phase reversals in the CORDIC pipeline.
4.  **Hardware Emulation**: High-fidelity Octave scripts that mirror RTL behavior for "Golden Reference" testing.