function writeLSSamples(n, matDir)

    outputDir = "samples/samplesLS_N"+n+".bin";
    disp(['Guardando: '+ outputDir]);
    load(matDir);

    nSamples = size(Ak,3);
    fileID = fopen(outputDir,'w');
    
    fwrite(fileID,n,'double');
    fwrite(fileID,nSamples,'double');
    fwrite(fileID,iterMINRES,'double');
    fwrite(fileID,tol,'double');

    for sample = 1:nSamples
        fwrite(fileID,reshape(Ak(:,:,sample)',1,[]),'double');
        fwrite(fileID,reshape(bk(:,:,sample)',1,[]),'double');
        fwrite(fileID,reshape(zk(:,:,sample),1,[]),'double');
    end
    fclose(fileID);
end



%{
function writeLSSamples(n, matDir)

    outputDir = "samples/samplesLS_N"+n+".txt";
    disp(['Guardando: '+ outputDir]);
    load(matDir);

    nSamples = size(Ak,3);
    fileID = fopen(outputDir,'w');
    fprintf(fileID,"N:\t%d\tnSamples:\t%d\n", n, nSamples);

    for sample = 1:nSamples
        fprintf(fileID,'%.30f\t',reshape(Ak(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(bk(:,:,sample)',1,[]));
        fprintf(fileID,'\n');
        fprintf(fileID,'%.30f\t',reshape(zk(:,:,sample),1,[]));
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end
%}