# High Level Synhesis, change prediction horizon
# unroll of matrix-vector multiplication external loop and vector operations
# ===============================================================================
# Alfonso Cortes Neira - Universidad Técnica Federico Santa María
# 01-10-2023
# Based on the work by Juan José Vásquez
# ===============================================================================

# Run Vitis HLS script with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls run_hls.tcl
# Execute Vitis HLS in shell mode with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls -i

set FORM "dense"
# Vector list
set N_HOR_LIST [list 8] 
# 2 3 4 4 5 6 7 8]
set UFACTOR_LIST [list 1 2 4 8 16 32 64]

set GIT_ROOT "C:/Users/Alfonso/Documents/GitHub/HLS_MPC"
set PRJ_ROOT "C:/dDesign/tesis/HLS/mpc"

# cd to script path
cd "C:/dDesign/tesis/HLS"
file mkdir "${PRJ_ROOT}/source"
file copy -force "${GIT_ROOT}/hls/source/mpc.cpp" "${PRJ_ROOT}/source/mpc.cpp"
file copy -force "${GIT_ROOT}/hls/source/mpc.hpp" "${PRJ_ROOT}/source/mpc.hpp"
file copy -force "${GIT_ROOT}/hls/source/utils.hpp" "${PRJ_ROOT}/source/utils.hpp"
file copy -force "${GIT_ROOT}/hls/source/system.hpp" "${PRJ_ROOT}/source/system.hpp"
file copy -force "${GIT_ROOT}/hls/source/MPCTestbench.cpp" "${PRJ_ROOT}/source/testbench.cpp"

foreach N_HOR $N_HOR_LIST {
    puts "\n////////////////////////////////"
	puts "\n//                            //"
	puts "\n// Creating N_HOR = ${N_HOR} project //"
	puts "\n//                            //"
	puts "\n////////////////////////////////\n"
	set PROBLEM "${FORM}_H${N_HOR}"
	# Create project.  -reset option allows to overwrite if the project exist
    open_project -reset mpc
    # Set HLS top function
    set_top mpc
	
    # Add sources, for the horizon
	file copy -force "${GIT_ROOT}/matlab/samples/MPC_motor_${FORM}_N${N_HOR}.cpp" "${PRJ_ROOT}/source/system.cpp"
    add_files "${PRJ_ROOT}/source/system.cpp" -cflags "-DHOR_SIZE=${N_HOR}" -csimflags "-DHOR_SIZE=${N_HOR}"
	add_files "${PRJ_ROOT}/source/system.hpp" -cflags "-DHOR_SIZE=${N_HOR}" -csimflags "-DHOR_SIZE=${N_HOR}"
	add_files "${PRJ_ROOT}/source/mpc.cpp" -cflags "-DHOR_SIZE=${N_HOR}" -csimflags "-DHOR_SIZE=${N_HOR}"
	add_files "${PRJ_ROOT}/source/mpc.hpp" -cflags "-DHOR_SIZE=${N_HOR}" -csimflags "-DHOR_SIZE=${N_HOR}"
	add_files "${PRJ_ROOT}/source/utils.hpp" -cflags "-DHOR_SIZE=${N_HOR}" -csimflags "-DHOR_SIZE=${N_HOR}"
	# Testbench source and Golden reference file. Note the -tb option
	file copy -force "${GIT_ROOT}/matlab/samples/MPC_motor_${FORM}_N${N_HOR}.bin" "${PRJ_ROOT}/source/samples.bin"
    add_files -tb "${PRJ_ROOT}/source/testbench.cpp" -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
    add_files -tb "${PRJ_ROOT}/source/samples.bin" -cflags " -Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
	
	foreach UFACTOR $UFACTOR_LIST {
		puts "\n--------------------------------\n"
		puts "\n   Solution Unroll_Factor = ${UFACTOR}   \n"
		puts "\n--------------------------------\n"
		set SOLUTION "multop_UN${UFACTOR}"
		# Set solution and flow target
		open_solution "${SOLUTION}" -flow_target vivado
		# Config solution with part for ZCU104 and 10ns target clock
		set_part {xczu7ev-ffvc1156-2-e}
		create_clock -period 10 -name default
		config_compile -unsafe_math_optimizations=true
		# Run C simulation
#		csim_design -clean
#		file copy -force "${PRJ_ROOT}/${SOLUTION}/csim/report/mpc_csim.log" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_csim.log"
		
		# Set pragmas with directives
		set_directive_array_partition -type complete -dim 2 "mvmult" A
		set_directive_array_partition -type complete -dim 1 "mvmult" B
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult" A
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult" R
		set_directive_unroll -factor ${UFACTOR} "mvmult/mvmult_row"
		set_directive_pipeline -II 1 "mvmult/mvmult_row"
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "vadd" A
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "vadd" B
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "vadd" R
		set_directive_unroll -factor ${UFACTOR} "vadd/vadd_row"
		set_directive_pipeline -II 1 "vadd/vadd_row"
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "vsub" A
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "vsub" B
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "vsub" R
		set_directive_unroll -factor ${UFACTOR} "vsub/vsub_row"
		set_directive_pipeline -II 1 "vsub/vsub_row"
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "max0" R
		set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "max0" A
		set_directive_unroll -factor ${UFACTOR} "max0/max0_row"
		set_directive_pipeline -II 1 "max0/max0_row"
		
		
		# Run synthesis
		csynth_design
		file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_csynth.txt"
	
		# Run cosimulation
		cosim_design
		file copy -force "${PRJ_ROOT}/${SOLUTION}/sim/report/verilog/mpc.log" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_sim.log"
		file copy -force "${PRJ_ROOT}/${SOLUTION}/sim/report/mpc_cosim.rpt" "${GIT_ROOT}/hls/results/${PROBLEM}_${SOLUTION}_sim.txt"
	
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