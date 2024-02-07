# High Level Synhesis, change number of constraints (DER)
# unroll of elemeto-to-element operations (add, sub, max)
# ===============================================================================
# Alfonso Cortes Neira - Universidad Técnica Federico Santa María
# 07-02-2024
# Based on the work by Juan José Vásquez
# ===============================================================================

# Run Vitis HLS script with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls der_op.tcl
# Execute Vitis HLS in shell mode with C:\Xilinx\Vitis_HLS\2022.2\bin\vitis_hls -i

set FORM "dense"
# Vector list
set CONST_LIST [list 8 12 20 32 16 24] 
# 2 3 4 5 6 7 8]
set UFACTOR_LIST [list 1 2 8 16]

set GIT_ROOT "C:/Users/Alfonso/Documents/GitHub/HLS_MPC"
set PRJ_ROOT "C:/dDesign/tesis/HLS/der"

# cd to script path
cd "C:/dDesign/tesis/HLS"
file mkdir "${PRJ_ROOT}/source"
file copy -force "${GIT_ROOT}/hls/source_der/mpc.cpp" "${PRJ_ROOT}/source/mpc.cpp"
file copy -force "${GIT_ROOT}/hls/source_der/mpc.hpp" "${PRJ_ROOT}/source/mpc.hpp"
file copy -force "${GIT_ROOT}/hls/source_der/utils.hpp" "${PRJ_ROOT}/source/utils.hpp"
file copy -force "${GIT_ROOT}/hls/source_der/system.hpp" "${PRJ_ROOT}/source/system.hpp"
file copy -force "${GIT_ROOT}/hls/source_der/MPCTestbench.cpp" "${PRJ_ROOT}/source/testbench.cpp"

foreach CONST $CONST_LIST {
    puts "\n////////////////////////////////"
	puts "\n//                            //"
	puts "\n// Creating CONST = ${CONST} project //"
	puts "\n//                            //"
	puts "\n////////////////////////////////\n"
	set PROBLEM "${FORM}_C${CONST}"
	# Create project.  -reset option allows to overwrite if the project exist
    open_project -reset der
    # Set HLS top function
    set_top mpc
	
    # Add sources, for the horizon
	file copy -force "${GIT_ROOT}/matlab/samples2/MPC_der_${FORM}_N2_C${CONST}.cpp" "${PRJ_ROOT}/source/system.cpp"
    add_files "${PRJ_ROOT}/source/system.cpp" -cflags "-DN_CONST=${CONST}" -csimflags "-DN_CONST=${CONST}"
	add_files "${PRJ_ROOT}/source/system.hpp" -cflags "-DN_CONST=${CONST}" -csimflags "-DN_CONST=${CONST}"
	add_files "${PRJ_ROOT}/source/mpc.cpp" -cflags "-DN_CONST=${CONST}" -csimflags "-DN_CONST=${CONST}"
	add_files "${PRJ_ROOT}/source/mpc.hpp" -cflags "-DN_CONST=${CONST}" -csimflags "-DN_CONST=${CONST}"
	add_files "${PRJ_ROOT}/source/utils.hpp" -cflags "-DN_CONST=${CONST}" -csimflags "-DN_CONST=${CONST}"
	# Testbench source and Golden reference file. Note the -tb option
	file copy -force "${GIT_ROOT}/matlab/samples2/MPC_der_${FORM}_N2_C${CONST}.bin" "${PRJ_ROOT}/source/samples.bin"
    add_files -tb "${PRJ_ROOT}/source/testbench.cpp" -cflags "-Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
    add_files -tb "${PRJ_ROOT}/source/samples.bin" -cflags " -Wno-unknown-pragmas" -csimflags "-Wno-unknown-pragmas"
	
	
	puts "\n--------------------------------\n"
	puts "\n   Solution sequential   \n"
	puts "\n--------------------------------\n"
	set SOLUTION "op_UN0"
	# Set solution and flow target
	open_solution "${SOLUTION}" -flow_target vivado
	# Config solution with part for ZCU104 and 10ns target clock
	set_part {xczu7ev-ffvc1156-2-e}
	create_clock -period 10 -name default
	config_compile -unsafe_math_optimizations=true
	# Run C simulation
#	csim_design -clean
#	file copy -force "${PRJ_ROOT}/${SOLUTION}/csim/report/mpc_csim.log" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_csim.log"
		
	# Set pragmas with directives
	set_directive_unroll -factor 1 "mvmult/mvmult_column"
	set_directive_pipeline -II 1 "mvmult/mvmult_column"
	set_directive_unroll -factor 1 "vadd/vadd_row"
	set_directive_pipeline -II 1 "vadd/vadd_row"
	set_directive_unroll -factor 1 "vsub/vsub_row"
	set_directive_pipeline -II 1 "vsub/vsub_row"
	set_directive_unroll -factor 1 "max0/max0_row"
	set_directive_pipeline -II 1 "max0/max0_row"

	# Run synthesis
	csynth_design
	file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_csynth.txt"
	
	# Run cosimulation
	cosim_design
	file copy -force "${PRJ_ROOT}/${SOLUTION}/sim/report/verilog/mpc.log" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_sim.log"
	file copy -force "${PRJ_ROOT}/${SOLUTION}/sim/report/mpc_cosim.rpt" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_sim.txt"
	
	# Set export settings. -version allows to avoid the Vitis HLS date error in export.
	config_export -display_name mpc -format ip_catalog -output "${GIT_ROOT}/hls/ips_der/der_${PROBLEM}_${SOLUTION}.zip" -rtl verilog -version 1.0
	# Run export + logic synthesis + implementation
	export_design -flow impl -rtl verilog -format ip_catalog -output "${GIT_ROOT}/hls/ips_der/der_${PROBLEM}_${SOLUTION}.zip"
	file copy -force "${PRJ_ROOT}/${SOLUTION}/impl/report/verilog/export_impl.rpt" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_impl.txt"
		
	# Move solution reports		
	file copy -force "${PRJ_ROOT}/${SOLUTION}/${SOLUTION}.log" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}.log"
		
	
	foreach UFACTOR $UFACTOR_LIST {
		puts "\n--------------------------------\n"
		puts "\n   Solution Unroll_Factor = ${UFACTOR}   \n"
		puts "\n--------------------------------\n"
		set SOLUTION "op_UN${UFACTOR}"
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
		set_directive_unroll "mvmult/mvmult_column"
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
		file copy -force "${PRJ_ROOT}/${SOLUTION}/syn/report/csynth.rpt" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_csynth.txt"
	
		# Run cosimulation
		cosim_design
		file copy -force "${PRJ_ROOT}/${SOLUTION}/sim/report/verilog/mpc.log" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_sim.log"
		file copy -force "${PRJ_ROOT}/${SOLUTION}/sim/report/mpc_cosim.rpt" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_sim.txt"
	
		# Set export settings. -version allows to avoid the Vitis HLS date error in export.
		config_export -display_name mpc -format ip_catalog -output "${GIT_ROOT}/hls/ips_der/der_${PROBLEM}_${SOLUTION}.zip" -rtl verilog -version 1.0
		# Run export + logic synthesis + implementation
		export_design -flow impl -rtl verilog -format ip_catalog -output "${GIT_ROOT}/hls/ips_der/der_${PROBLEM}_${SOLUTION}.zip"
		file copy -force "${PRJ_ROOT}/${SOLUTION}/impl/report/verilog/export_impl.rpt" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}_impl.txt"
		
		# Move solution reports		
		file copy -force "${PRJ_ROOT}/${SOLUTION}/${SOLUTION}.log" "${GIT_ROOT}/hls/results_der/${PROBLEM}_${SOLUTION}.log"
		
	}
}

puts "\nDONE\n"
set exit [gets stdin]