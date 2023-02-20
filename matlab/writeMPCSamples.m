function writeMPCSamples(n, matDir)

    outputDir = "samples/samplesMPC_N"+n+".bin";

    disp(['Guardando: '+ outputDir]);
    load(matDir);

    I = size(M_hat,1);
    nSamples = size(u_star,2);

    fileID = fopen(outputDir,'w');
    fwrite(fileID,n,'double');
    fwrite(fileID,I,'double');
    fwrite(fileID,nSamples,'double');
    fwrite(fileID,iterPDIP,'double');
    fwrite(fileID,iterMINRES,'double');
    fwrite(fileID,tol,'double');
    

    fwrite(fileID,reshape(xmin',1,[]),'double');
    fwrite(fileID,reshape(xmax',1,[]),'double');
    fwrite(fileID,reshape(umin',1,[]),'double');
    fwrite(fileID,reshape(umax',1,[]),'double');
    fwrite(fileID,reshape(Acal',1,[]),'double');
    fwrite(fileID,reshape(AcalQOcal',1,[]),'double');
    fwrite(fileID,reshape(H',1,[]),'double');
    fwrite(fileID,reshape(M_hat',1,[]),'double');
    fwrite(fileID,reshape(L_invLast',1,[]),'double');
    fwrite(fileID,reshape(A',1,[]),'double');
    fwrite(fileID,reshape(B',1,[]),'double');
     
    for sample = 1:nSamples
        fwrite(fileID,reshape(x(:,sample)',1,[]),'double');
        fwrite(fileID,reshape(r(:,sample)',1,[]),'double');
        fwrite(fileID,reshape(u_star(:,sample),1,[]),'double');
    end
    fclose(fileID);
end

