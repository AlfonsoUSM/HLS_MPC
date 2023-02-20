function writePDIPSamples(n, matDir)
    outputDir = "samples/samplesPDIP_N"+n+".bin";
    disp(['Guardando: '+ outputDir]);
    load(matDir);

    I = size(M_hat,1);
    nSamples = size(u_tilde_star,3);

    fileID = fopen(outputDir,'w');
    fwrite(fileID,n,'double');
    fwrite(fileID,I,'double');
    fwrite(fileID,nSamples,'double');
    fwrite(fileID,iterPDIP,'double');
    fwrite(fileID,iterMINRES,'double');
    fwrite(fileID,tol,'double');
    fwrite(fileID,reshape(H',1,[]),'double');
    fwrite(fileID,reshape(M_hat',1,[]),'double');

    for sample = 1:nSamples
        fwrite(fileID,reshape(h_tilde(:,:,sample)',1,[]),'double');
        fwrite(fileID,reshape(c_hat(:,:,sample)',1,[]),'double');
        fwrite(fileID,reshape(u_tilde_star(:,:,sample),1,[]),'double');
    end
    fclose(fileID);
end


%{
function writePDIPSamples(n, matDir)
    outputDir = "samples/samplesPDIP_N"+n+".txt";
    disp(['Guardando: '+ outputDir]);
    load(matDir);

    I = size(Mx,1);
    nSamples = size(u_tilde,3);

    fileID = fopen(outputDir,'w');
    fprintf(fileID,"N:\t%d\tI:\t%d\tnSamples:\t%d\n", n, I, nSamples);
    fprintf(fileID,'%.30f\t',reshape(H',1,[]));
    fprintf(fileID,'\n');    
    fprintf(fileID,'%.30f\t',reshape(Mx',1,[]));
    fprintf(fileID,'\n'); 

    for sample = 1:nSamples
        fprintf(fileID,'%.30f\t',reshape(h(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(cx(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(u_tilde(:,:,sample),1,[]));
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end
%}