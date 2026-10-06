# ==============================================================================
# Vivado Project Generation Script (Project Mode)
# Script Location: scripts/range_detector.tcl
# ==============================================================================

# Determine paths relative to the scripts/ folder
set script_dir [file dirname [file normalize [info script]]]
set module_dir [file normalize "$script_dir/.."]

# Project Configuration
set proj_name        "range_detector_proj"
set target_part      "xc7z020clg400-1"
set top_design_name  "range_detector"
set top_tb_name      "tb_range_detector"
set proj_dir         "$module_dir/project"

# 1. Create Vivado Project
create_project $proj_name $proj_dir -part $target_part -force

# 2. Configure Language Properties
set_property target_language VHDL [current_project]
set_property simulator_language VHDL [current_project]

# 3. Import Explicit IP Cores
set target_ip_files [list \
    "$module_dir/ipcores/adc_fifo/adc_fifo.xci" \
    "$module_dir/ipcores/clk_dsp/clk_dsp.xci" \
]

puts "--> Adding targeted IP cores..."
foreach ip_file $target_ip_files {
    if {[file exists $ip_file]} {
        read_ip $ip_file
        generate_target all [get_files $ip_file]
        puts "    Added IP: $ip_file"
    } else {
        puts "    \[WARNING\] IP file not found: $ip_file"
    }
}

# 4. Import Explicit Design Sources (rtl_src)
set target_rtl_files [list \
    "$module_dir/rtl_src/radar_pkg.vhd" \
    "$module_dir/rtl_src/DSP_wrapper.vhd" \
    "$module_dir/rtl_src/fir_impl.vhd" \
    "$module_dir/rtl_src/complex_convolution.vhd" \
    "$module_dir/rtl_src/range_detector.vhd" \
]

puts "--> Adding targeted HDL sources..."
foreach rtl_file $target_rtl_files {
    if {[file exists $rtl_file]} {
        add_files -fileset sources_1 $rtl_file
        puts "    Added RTL: $rtl_file"
    } else {
        puts "    \[WARNING\] RTL file not found: $rtl_file"
    }
}

# 5. Import Explicit Simulation Sources & Waveform Config (rtl_sim)
set target_sim_files [list \
    "$module_dir/rtl_sim/rtl_golden_ref_vector.vhd" \
    "$module_dir/rtl_sim/tb_range_detector.vhd" \
    "$module_dir/rtl_sim/range_detector.wcfg" \
]

puts "--> Adding targeted simulation sources..."
foreach sim_file $target_sim_files {
    if {[file exists $sim_file]} {
        add_files -fileset sim_1 $sim_file
        puts "    Added Sim Asset: $sim_file"
    } else {
        puts "    \[WARNING\] File not found: $sim_file"
    }
}

# 6. Import Constraint Files (.xdc)
if {[file exists "$module_dir/constraints"]} {
    puts "--> Adding constraints from $module_dir/constraints..."
    add_files -fileset constrs_1 "$module_dir/constraints"
}

# 7. Update Compile Order and Set Top Modules
update_compile_order -fileset sources_1
update_compile_order -fileset sim_1

if {[get_filesets sources_1] ne ""} {
    set_property top $top_design_name [get_filesets sources_1]
    update_compile_order -fileset sources_1
}

if {[get_filesets sim_1] ne ""} {
    set_property top $top_tb_name [get_filesets sim_1]
    update_compile_order -fileset sim_1
}

puts "=========================================================================="
puts " Project '$proj_name' created successfully at:"
puts " $proj_dir"
puts " Top Level Entity: $top_design_name"
puts " Simulation Top:   $top_tb_name"
puts "=========================================================================="