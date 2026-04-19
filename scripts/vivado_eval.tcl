if { $argc < 6 } {
  puts "Usage: vivado -mode batch -source scripts/vivado_eval.tcl -tclargs <verilog> <top> <part> <flow_mode> <period_ns> <report_dir>"
  exit 1
}

set verilog_file [lindex $argv 0]
set top_module [lindex $argv 1]
set part_name [lindex $argv 2]
set flow_mode [lindex $argv 3]
set period_ns [lindex $argv 4]
set report_dir [lindex $argv 5]

file mkdir $report_dir

set util_rpt [file join $report_dir "utilization.rpt"]
set timing_rpt [file join $report_dir "timing.rpt"]
set power_rpt [file join $report_dir "power.rpt"]

read_verilog $verilog_file
synth_design -top $top_module -part $part_name -mode out_of_context

# Deterministic combinational timing requirement.
set_max_delay $period_ns -from [all_inputs] -to [all_outputs]

if { $flow_mode eq "post_route" } {
  opt_design
  place_design
  route_design
} elseif { $flow_mode ne "synth_only" } {
  puts "Unsupported flow_mode '$flow_mode' (expected synth_only or post_route)"
  exit 2
}

report_utilization -file $util_rpt
report_timing_summary -delay_type max -max_paths 1 -file $timing_rpt
report_power -file $power_rpt

exit 0
