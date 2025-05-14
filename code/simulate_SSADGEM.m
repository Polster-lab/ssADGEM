%load('Human-GEM/model/Human-GEM.mat')
load('../models/ROSMAP_GEM.mat')
temp = setHamsMedia(model);
[solutions, goodRxns] = randomSampling(temp,10000,false,true,true,[],true);
goodRxns = model.rxns(goodRxns);
annotation  = readtable('../data/ROSMAP_annotation_processed.txt','Delimiter','\t');
ADpos = find(strcmpi(annotation.AD,'AD'));
annotation(138,:) = [];
annotation(140,:) = [];

b = height(annotation);
baseModel = ihuman;

for i=1:1
    %data = readtable(['../results/DE_analysis/DE_genes_RM_' num2str(i) '.txt'],'Delimiter','\t');
    patient = annotation.patient{i};
    disp(['cluster #' num2str(i)])

    %try
        load(['../models/ROSMAP_AD_cluster_' num2str(i) '.mat'],'clusterModel')
        temp = setHamsMedia(clusterModel);
        %identify good rxns for random sampling
        [ia,ib] = ismember(goodRxns,temp.rxns);
        solutions = randomSampling(temp,10000,false,true,true,ib(ib>0),true);
            
    %catch
        disp(['Model for cluster: ' num2str(i) ' not found'])
    %end
end
haveFlux = abs(solutions)>1E-3;
activeRxns = (sum(haveFlux,2));
