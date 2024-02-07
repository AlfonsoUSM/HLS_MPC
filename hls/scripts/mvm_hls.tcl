# High Level Synhesis, matrix-vector multiplication
# change matrix size and unroll factor, row-wise and column-wise
# ===============================================================================
# Alfonso Cortes Neira - Universidad Técnica Federico Santa María
# 13-01-2024
# ===============================================================================

# Run Vitis HLS script with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls hls_mvm.tcl
# Execute Vitis HLS in shell mode with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls -i

set MRATIO_LIST [list 1 2 4 8 16 32]
set UFACTOR_LIST [list 1 2 4 8 16 32 64]

set GIT_ROOT "C:/Users/Alfonso/Documents/GitHub/HLS_MPC"
set PRJ_ROOT "C:/dDesign/tesis/HLS/mvm"

# cd to script path
cd "C:/dDesign/tesis/HLS"
file mkdir "${PRJ_ROOT}/source"
file copy -force "${GIT_ROOT}/hls/source_mvm/mvmult.cpp" "${PRJ_ROOT}/source/mvmult.cpp"
file copy -force "${GIT_ROOT}/hls/source_mvm/mvmult.hpp" "${PRJ_ROOT}/source/mvmult.hpp"


foreach RATIO $MRATIO_LIST {
    puts "\n////////////////////////////////"
	puts "\n//                            //"
	puts "\n// Creating RATIO = ${RATIO} project //"
	puts "\n//                            //"
	puts "\n////////////////////////////////\n"
	set N_SIZE [expr {2 * $RATIO}]
	set M_SIZE [expr {64 / $RATIO}]
	# Create project.  -reset option allows to overwrite if the project exist
    open_project -reset mvm
	
    # Add sources
    add_files "${PRJ_ROOT}/source/mvmult.cpp" -cflags "-DPARAM_N=${N_SIZE} -DPARAM_M=${M_SIZE}" -csimflags "-DPARAM_N=${N_SIZE} -DPARAM_M=${M_SIZE}"
	add_files "${PRJ_ROOT}/source/mvmult.hpp" -cflags "-DPARAM_N=${N_SIZE} -DPARAM_M=${M_SIZE}" -csimflags "-DPARAM_N=${N_SIZE} -DPARAM_M=${M_SIZE}"
	
    # Set HLS top function column wise
    set_top mvmult_rowwise
	set PROBLEM "${N_SIZE}x${M_SIZE}_rowwise_in2"
	puts "\n--------------------------------\n"
	puts "\n    Problem = ${PROBLEM}    \n"
	
	foreach UFACTOR $UFACTOR_LIST {
		if { $UFACTOR < $M_SIZE } {
			puts "\n--------------------------------\n"
			puts "\n   Solution Unroll_Factor = ${UFACTOR}   \n"
			puts "\n--------------------------------\n"
			set SOLUTION "UN${UFACTOR}"
			# Set solution and flow target
			open_solution "${SOLUTION}" -flow_target vivado
			# Config solution with part for ZCU104 and 10ns target clock
			set_part {xczu7ev-ffvc1156-2-e}
			create_clock -period 10 -name default
			config_compile -unsafe_math_optimizations=true
		
			# Set pragmas with directives
			set_directive_array_partition -dim 2 -factor ${UFACTOR} -type cyclic "mvmult_rowwise" A
			set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult_rowwise" B
#			set_directive_loop_flatten -off "mvmult_rowwise/outer_loop"
			set_directive_unroll -factor ${UFACTOR} "mvmult_rowwise/inner_loop"
			set_directive_pipeline -II 1 "mvmult_rowwise/inner_loop"
		
			# Run synthesis
			csynth_design
			file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}_csynth.txt"
			# Move solution reports		
#			file copy -force "${PRJ_ROOT}/${SOLUTION}/${SOLUTION}.log" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}.log"
		}
	}
	
#	set PROBLEM "${N_SIZE}x${M_SIZE}_rowwise_out"
#	puts "\n--------------------------------\n"
#	puts "\n    Problem = ${PROBLEM}    \n"
#	
#	foreach UFACTOR $UFACTOR_LIST {
#		if { $UFACTOR <= $N_SIZE } {
#			puts "\n--------------------------------\n"
#			puts "\n   Solution Unroll_Factor = ${UFACTOR}   \n"
#			puts "\n--------------------------------\n"
#			set SOLUTION "UN${UFACTOR}"
#			# Set solution and flow target
#			open_solution "${SOLUTION}" -flow_target vivado
#			# Config solution with part for ZCU104 and 10ns target clock
#			set_part {xczu7ev-ffvc1156-2-e}
#			create_clock -period 10 -name default
#			config_compile -unsafe_math_optimizations=true
#		
#			# Set pragmas with directives
#			set_directive_array_partition -type complete -dim 2 "mvmult_rowwise" A
#			set_directive_array_partition -type complete -dim 1 "mvmult_rowwise" B
#			set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult_rowwise" A
#			set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult_rowwise" R
#			set_directive_unroll -factor ${UFACTOR} "mvmult_rowwise/outer_loop"
#			set_directive_pipeline -II 1 "mvmult_rowwise/outer_loop"
#		
#			# Run synthesis
#			csynth_design
#			file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}_csynth.txt"
#			# Move solution reports		
#			file copy -force "${PRJ_ROOT}/${SOLUTION}/${SOLUTION}.log" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}.log"
#		}
#	}

    set_top mvmult_columnwise
	set PROBLEM "${N_SIZE}x${M_SIZE}_columnwise_in2"
	puts "\n--------------------------------\n"
	puts "\n    Problem = ${PROBLEM}    \n"
	
	foreach UFACTOR $UFACTOR_LIST {
		if { $UFACTOR < $N_SIZE } {
			puts "\n--------------------------------\n"
			puts "\n   Solution Unroll_Factor = ${UFACTOR}   \n"
			puts "\n--------------------------------\n"
			set SOLUTION "UN${UFACTOR}"
			# Set solution and flow target
			open_solution "${SOLUTION}" -flow_target vivado
			# Config solution with part for ZCU104 and 10ns target clock
			set_part {xczu7ev-ffvc1156-2-e}
			create_clock -period 10 -name default
			config_compile -unsafe_math_optimizations=true
		
			# Set pragmas with directives
			set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult_columnwise" A
			set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult_columnwise" R
#			set_directive_loop_flatten -off "mvmult_columnwise/outer_loop"
			set_directive_unroll -factor ${UFACTOR} "mvmult_columnwise/inner_loop"
			set_directive_pipeline -II 1 "mvmult_columnwise/inner_loop"
		
			# Run synthesis
			csynth_design
			file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}_csynth.txt"
			# Move solution reports		
#			file copy -force "${PRJ_ROOT}/${SOLUTION}/${SOLUTION}.log" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}.log"
		}
	}
	
#	set PROBLEM "${N_SIZE}x${M_SIZE}_columnwise_out"
#	puts "\n--------------------------------\n"
#	puts "\n    Problem = ${PROBLEM}    \n"
#	
#	foreach UFACTOR $UFACTOR_LIST {
#		if { $UFACTOR <= $M_SIZE } {
#			puts "\n--------------------------------\n"
#			puts "\n   Solution Unroll_Factor = ${UFACTOR}   \n"
#			puts "\n--------------------------------\n"
#			set SOLUTION "UN${UFACTOR}"
#			# Set solution and flow target
#			open_solution "${SOLUTION}" -flow_target vivado
#			# Config solution with part for ZCU104 and 10ns target clock
#			set_part {xczu7ev-ffvc1156-2-e}
#			create_clock -period 10 -name default
#			config_compile -unsafe_math_optimizations=true
#		
#			# Set pragmas with directives
#			set_directive_array_partition -type complete -dim 1 "mvmult_columnwise" A
#			set_directive_array_partition -type complete -dim 1 "mvmult_columnwise" R
#			set_directive_array_partition -dim 2 -factor ${UFACTOR} -type cyclic "mvmult_columnwise" A
#			set_directive_array_partition -dim 1 -factor ${UFACTOR} -type cyclic "mvmult_columnwise" B
#			set_directive_unroll -factor ${UFACTOR} "mvmult_columnwise/outer_loop"
#			set_directive_pipeline -II 1 "mvmult_columnwise/outer_loop"
#		
#			# Run synthesis
#			csynth_design
#			file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}_csynth.txt"
#			# Move solution reports		
#			file copy -force "${PRJ_ROOT}/${SOLUTION}/${SOLUTION}.log" "${GIT_ROOT}/hls/results_mvm/${PROBLEM}_${SOLUTION}.log"
#		}
#	}	
}

puts "\n    DONE    \n"
set exit [gets stdin]