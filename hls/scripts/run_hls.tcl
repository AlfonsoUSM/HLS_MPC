# High Level Synhesis, change prediction horizon
# ===============================================================================
# Alfonso Cortes Neira - Universidad Técnica Federico Santa María
# 12-09-2023
# Based on the work by Juan José Vásquez
# ===============================================================================

# Run Vitis HLS script with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls run_hls.tcl
# Execute Vitis HLS in shell mode with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls -i

set FORM "dense"
# Vector list
set N_HOR_LIST [list 4] 
# 4 8 16 32 64]
set UFACTOR_LIST [list 1]

set GIT_ROOT "C:/Users/Alfonso/Documents/GitHub/HLS_MPC"
set PRJ_ROOT "C:/dDesign/tesis/HLS/mpc"

# cd to script path
cd "C:/dDesign/tesis/HLS"
file mkdir "${PRJ_ROOT}/source"
file copy -force "${GIT_ROOT}/hls/source/mpc.cpp" "${PRJ_ROOT}/source/mpc.cpp"
file copy -force "${GIT_ROOT}/hls/source/mpc.hpp" "${PRJ_ROOT}/source/mpc.hpp"
file copy -force "${GIT_ROOT}/hls/source/utils.hpp" "${PRJ_ROOT}/source/utils.hpp"
file copy -force "${GIT_ROOT}/hls/source/system.hpp" "${PRJ_ROOT}/source/system.hpp"

foreach N_HOR $N_HOR_LIST {

    puts "\n\n\nCreating N_HOR = ${N_HOR} project \n\n\n"
	set PROBLEM "${FORM}_H${N_HOR}"
	# Create project.  -reset option allows to overwrite if the project exist
    open_project -reset mpc
    # Set HLS top function
    set_top mpc
	
	# Modify N_HOR define and sample
	
    # Add config file to sources
    add_files "${PRJ_ROOT}/source/mpc.cpp ${PRJ_ROOT}/source/mpc.hpp ${PRJ_ROOT}/source/utils.hpp ${PRJ_ROOT}/source/system.hpp" 
	file copy -force "${GIT_ROOT}/matlab/samples/MPC_motor_${FORM}_N${N_HOR}.cpp" "${PRJ_ROOT}/source/system.cpp"
    add_files "${PRJ_ROOT}/source/system.cpp"
	# Testbench source. Note the -tb option
    add_files -tb "${GIT_ROOT}/hls/source/MPCTestbench.cpp" -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
    # Golden reference file
    add_files -tb "${GIT_ROOT}/matlab/samples/MPC_motor_${FORM}_N${N_HOR}.bin" -cflags " -Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
	
	foreach UFACTOR $UFACTOR_LIST {
		set SOLUTION "op_UN${UFACTOR}"
		# Set solution and flow target
		open_solution "${SOLUTION}" -flow_target vivado
		# Run C simulation
		#csim_design -clean
		#file copy -force "${PRJ_ROOT}/${SOLUTION}/csim/report/mpc_csim.log" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_csim.log"
		
		# Config solution with part for ZCU104 and 10ns target clock
		set_part {xczu7ev-ffvc1156-2-e}
		create_clock -period 10 -name default
		# Run synthesis
		csynth_design
		file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_csynth.txt"
	
		# Run cosimulation
		cosim_design
		file copy -force "${PRJ_ROOT}/${SOLUTION}/sim/report/verilog/mpc.log" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_sim.log"
	
		# Set export settings. -version allows to avoid the Vitis HLS date error in export.
		config_export -display_name mpc -format ip_catalog -output "${GIT_ROOT}/hls/ips/${PROBLEM}_${SOLUTION}.zip" -rtl verilog -version 1.0
		# Run export + logic synthesis + implementation
		export_design -flow impl -rtl verilog -format ip_catalog -output "${GIT_ROOT}/hls/ips/${PROBLEM}_${SOLUTION}.zip"
		file copy -force "${PRJ_ROOT}/${SOLUTION}/impl/report/verilog/export_impl.rpt" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_impl.txt"
		
		# Move solution reports		
		file copy -force "${PRJ_ROOT}/${SOLUTION}/${SOLUTION}.log" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}.log"
		
	}
}

set exit [gets stdin]