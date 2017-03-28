nPerAxisList=4*floor(2.^[1:0.5:8]);
for nPerAxis=nPerAxisList
    solveAndPlot_Lshape_pde_shishkin_batch;
    save(['Lshape_pde_shishkin3_',num2str(nPerAxis)],'basis','fun2id','id2fun','mesh0','meshType','n','w','nPerAxis'...
                                                   ,'numSol','S','u','vecf');
    
    solveAndPlot_Lshape_pde_uniform_batch;
    save(['Lshape_pde_uniform2_',num2str(nPerAxis)],'basis','fun2id','id2fun','mesh0','meshType','n','w','nPerAxis'...
                                                   ,'numSol','S','u','vecf');
end