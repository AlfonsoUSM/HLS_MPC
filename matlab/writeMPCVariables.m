%% ADMM MPC Global Variables
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 25-07-2023
% Write variables with float literals 
% ===============================================================================

function writeMPCVariables(n, matDir, plant, data_t)
    cpp_file = 'test.cpp';
    fID = fopen(cpp_file, 'w');
    %load(matDir);

% H
    H_qp = full(H);
    fprintf(fID, 'data_t H_qp[N_QP][N_QP] =  {\n');
    for i = 1:N_QP
        for j = 1:N_QP
            fprintf(fID, cpp_float(H_qp(i,j)));
            fprintf(fID, ',\t');
        end
        fprintf(fID, '\n');
    end
    fprintf(fID, '};\n\n');

% R_inv
    fprintf(fID, 'data_t R_inv[N_QP][N_QP] =  {\n');
    for i = 1:N_QP
        for j = 1:N_QP
            fprintf(fID, cpp_float(R_inv(i,j)));
            fprintf(fID, ',\t');
        end
        fprintf(fID, '\n');
    end
    fprintf(fID, '};\n\n');

    fclose(fID);
end