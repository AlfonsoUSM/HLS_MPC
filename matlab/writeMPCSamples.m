function writeMPCSamples(n, matDir, plant, data_t)

    binfile = "samples/MPC_"+plant+"_N"+n+"_"+data_t+".bin";
%     txtfile = "samples/samplesMPC_N"+n+".txt";

    disp('Guardando: '+ binfile);
    load(matDir);

    nSamples = size(uk,2);
    H_qp = full(H);
    C_qp = full(M_hat);
    RMN = full(RhoMt_neg);

    binfileID = fopen(binfile,'w');


%     txtfileID = fopen(txtfile,'w');
%     % H
%     fprintf(txtfileID, 'data_t H_qp[N_QP][N_QP] =\t{');
%     for row = 1:(N_QP-1)
%         fprintf(txtfileID, '{');
%         fprintf(txtfileID, '%f, ', H_qp(row,1:(N_QP-1)));
%         fprintf(txtfileID, '%f},\n\t\t\t\t', H_qp(row,N_QP));
%     end
%     fprintf(txtfileID, '{');
%     fprintf(txtfileID, '%f, ', H_qp(N_QP,1:(N_QP-1)));
%     fprintf(txtfileID, '%f}\n\n', H_qp(N_QP,N_QP));
%     % h
%     fprintf(txtfileID, 'data_t h_qp[N_QP] =\t{');
%     fprintf(txtfileID, '%f, ', h(1:(N_QP-1)));
%     fprintf(txtfileID, '%f}\n\n', H_qp(row,N_QP));
%     % M_qp
%     fprintf(txtfileID, 'data_t C_qp[N_QP][N_QP] =\t{');
%     for row = 1:(M_QP-1)
%         fprintf(txtfileID, '{');
%         fprintf(txtfileID, '%f, ', C_qp(row,1:(N_QP-1)));
%         fprintf(txtfileID, '%f},\n\t\t\t\t', C_qp(row,N_QP));
%     end
%     fprintf(txtfileID, '{');
%     fprintf(txtfileID, '%f, ', C_qp(M_QP,1:(N_QP-1)));
%     fprintf(txtfileID, '%f}\n\n', C_qp(M_QP,N_QP));
%     % g
%     fprintf(txtfileID, 'data_t g[(2*N_QP)] =\t{');
%     fprintf(txtfileID, '%f, ', g(1:(2*N_QP-1)));
%     fprintf(txtfileID, '%f}\n\n', g(row,2*N_QP));
%     fclose(txtfileID);

    fwrite(binfileID, n,'uint8');
    fwrite(binfileID, N_SYS,'uint8');
    fwrite(binfileID, M_SYS,'uint8');
    fwrite(binfileID, P_SYS,'uint8');
    fwrite(binfileID, nSamples,'uint64');
    fwrite(binfileID, N_QP,'uint8');
    fwrite(binfileID, M_QP,'uint8');
    fwrite(binfileID, iter,'uint16');  

    fwrite(binfileID,reshape(xmin',1,[]),data_t);
    fwrite(binfileID,reshape(xmax',1,[]),data_t);
    fwrite(binfileID,reshape(umin',1,[]),data_t);
    fwrite(binfileID,reshape(umax',1,[]),data_t);
    fwrite(binfileID,reshape(H_qp',1,[]),data_t);
    fwrite(binfileID,reshape(h',1,[]),data_t);
    fwrite(binfileID,reshape(C_qp',1,[]),data_t);
    fwrite(binfileID,reshape(g',1,[]),data_t);
    fwrite(binfileID,reshape(A',1,[]),data_t);
    fwrite(binfileID,reshape(B',1,[]),data_t);
    fwrite(binfileID,Rho,data_t);
    fwrite(binfileID,reshape(R_inv',1,[]),data_t);
    fwrite(binfileID,reshape(RMN',1,[]),data_t);
     
    for sample = 1:nSamples
%         fwrite(fileID,reshape(c_hat(:,sample)',1,[]),'double');
        fwrite(binfileID,reshape(xk(:,sample)',1,[]),data_t);
        fwrite(binfileID,reshape(uk(:,sample)',1,[]),data_t);
%         fwrite(fileID,reshape(theta(:,sample),1,[]),'double');
    end

    fclose(binfileID);
%     fclose(txtfileID);
end

