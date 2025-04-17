%Get latest version of humanGEM and load xml file with 
git clone --depth=1 https://github.com/SysBioChalmers/Human-GEM.git
ihuman = load("Human-GEM/model/Human-GEM.mat");
ihuman = ihuman.ihuman;
mkdir('../results')
modelVer = ihuman.version;
dateStr  = datetime('today');
SSmodels = false;
%%Generate pan-ROSMAP model
geneCounts = readtable('../data/ROSMAP_meanCounts_processed_all.txt','Delimiter','\t');
[~,genes2keep] = ismember(ihuman.genes,table2cell(geneCounts(:,1))); 
geneCounts = geneCounts(genes2keep(genes2keep>0),:);
data_struct = struct();
%From R script: #normalized to CPM and keep those genes with  >=1 CPM in at least one sample (permissive approach)
values = table2array(geneCounts(:,2));
idxs   = find(values>0);
%set data structure for ftINIT
data_struct.genes = table2cell(geneCounts(idxs,1));
data_struct.levels = values(idxs,:);
data_struct.tissues = {'ROSMAP'};
data_struct.threshold  = 1;
cd Human-GEM/code/tINIT/
%prepDataHumanGEM = prepHumanModelForftINIT(ihuman,false,'../../data/metabolicTasks/metabolicTasks_Essential.txt','../../model/reactions.tsv');
[model,~,addedRxnsForTasks_ROSMAP, deletedRxnsInINIT_ROSMAP] = ftINIT(prepDataHumanGEM, data_struct.tissues{1}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), false, true);
model.id = 'Human-GEM pan-ROSMAP';
cd ../../../
save('../models/ROSMAP_GEM.mat','model')
writeYAMLmodel(model,'../models/ROSMAP_GEM.yml')
ROSMAP = model;
%%Generate NCI model
allCounts   = readtable('../data/ROSMAP_counts_processed_all.txt','Delimiter','\t');
annotation  = readtable('../data/ROSMAP_annotation_processed.txt','Delimiter','\t');
allCounts = allCounts(genes2keep(genes2keep>0),:);
values      = table2array(allCounts(:,1:(end-1)));
[a,nSamples] = size(values);
normVals    = values;
newVals = table2array(geneCounts(:,2));
%normalize counts (cpm)
for i=1:nSamples
    normVals(:,i) = values(:,i)./(sum(values(:,i),1)/1E6);
    normVals(normVals(:,i)<1,i) = 0;
end
baseModel = prepDataHumanGEM.refModel;
% now build a matrix to compare presence/absence of reactions
NCI     = find(strcmp(annotation.AD,'No_AD'));
NCIcounts = normVals(:,NCI);
logiNCI = logical(NCIcounts);
sumas = sum(logiNCI,2);
%Keep genes expressed in at least 50% of NCI samples
idxsSum = find(sumas>=0.5*length(NCI));
tissues = {'NCI'};
data_struct.tissues = tissues;
countsMatrix = normVals(idxsSum);
data_struct.genes = geneCounts.gene(idxsSum);
data_struct.levels  = [newVals(idxsSum)];
data_struct.threshold  = 1;
cd Human-GEM/code/tINIT/
try
    [NCI_model, ~, addedRxnsForTasks_NCI, deletedRxnsInINIT_NCI]= ftINIT(prepDataHumanGEM,'NCI', [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
    NCI_model.id = 'Human-GEM NCI-ROSMAP';
    cd ../../../
    save('../models/ROSMAP_NCI_GEM.mat','NCI_model')
    writeYAMLmodel(model,'../models/ROSMAP_GEM_NCI.yml')
catch
end
%%Get sample-specific models
if SSmodels
    ADpos = find(strcmp(annotation.AD,'AD'));
    cd Human-GEM/code/tINIT/
    patients = allCounts.Properties.VariableNames;
    patients = strrep(patients,'x','');
    %prepDataHumanGEM2 = prepHumanModelForftINIT(model,false,'../../data/metabolicTasks/metabolicTasks_Essential.txt','../../model/reactions.tsv');
    %pseudosample = mean(normVals(:,NCI),2);
    compRxnMat= zeros(length(baseModel.rxns), nSamples);
    compgnsMat= zeros(length(baseModel.genes), nSamples);
    compMetMat= zeros(length(baseModel.mets), nSamples);
    compTaskAdd= zeros(length(baseModel.rxns), nSamples);
    compINITrmv= zeros(length(baseModel.rxns), nSamples);
    %parpool(8)
    nonExpressed = cell(1,numel(genes2keep));
    addedRxns_SS = cell(nSamples,1);
    deletedRxns_SS = cell(nSamples,1);
    for i=1:nSamples
        %cd Human-GEM/code/tINIT/
        patient = annotation.patient{i};
        disp([num2str(i) ': ' patient])
        data_struct = struct();
        data_struct.genes = geneCounts.gene;
        index = find(strcmp(patients,patient));
        data_struct.levels  = normVals;
        data_struct.tissues = annotation.patient;
        data_struct.threshold  = 1;
        try
            [sampleModel, ~, addedRxnsForTasks, deletedRxnsInINIT]= ftINIT(prepDataHumanGEM, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
            save(['../../../../models/' patient '.mat'],'sampleModel')
            %    load(['../models_pseudo/' patient '.mat'],'sampleModel')
            sampleModel.id = ['ROSMAP_' patient];
            [~,rxnPresence] = ismember(sampleModel.rxns,baseModel.rxns);
            compRxnMat(rxnPresence,i) = 1;
            [~,gnsPresence] = ismember(sampleModel.genes,baseModel.genes);
            compgnsMat(gnsPresence,i) = 1;
            [~,metPresence] = ismember(sampleModel.mets,baseModel.mets);
            compMetMat(metPresence,i) = 1;
            [~,gnsPresence] = ismember(sampleModel.genes,baseModel.genes);
            compgnsMat(gnsPresence,i) = 1;
            [~,rxnPresence] = ismember(addedRxnsForTasks,baseModel.rxns);
            compTaskAdd(rxnPresence,i) = 1;
            [~,rxnPresence] = ismember(deletedRxnsInINIT,baseModel.rxns);
            compINITrmv(rxnPresence,i) = 1;
            addedRxns_SS(i,1) = {addedRxnsForTasks};
            deletedRxns_SS(i,1) = {deletedRxnsInINIT};
        catch
        end
        clc
    end
end
%%Get AD cluster models
%get GRN cluster info
cluster1 = readtable(['../data/samples_GRN_cluster_1.txt'],'Delimiter',',','NumHeaderLines',0);
cluster1 = strrep(cluster1.Properties.VariableNames,'x','');
cluster2 = readtable(['../data/samples_GRN_cluster_2.txt'],'Delimiter',',','NumHeaderLines',0);
cluster2 = strrep(cluster2.Properties.VariableNames,'x','');
cluster3 = readtable(['../data/samples_GRN_cluster_3.txt'],'Delimiter',',','NumHeaderLines',0);
cluster3 = strrep(cluster3.Properties.VariableNames,'x','');
annotation.cluster = zeros(height(annotation),1);
[ia,ib] = ismember(cluster1,annotation.patient);
annotation.cluster(ib(ib>0)) = 1;
[ia,ib] = ismember(cluster2,annotation.patient);
annotation.cluster(ib(ib>0)) = 2;
[ia,ib] = ismember(cluster3,annotation.patient);
annotation.cluster(ib(ib>0)) = 3;
cd Human-GEM/code/tINIT/
for i=1:3
    patient_idxs =find(annotation.cluster == i);
    disp(['cluster' num2str(i)])
    data_struct = struct();
    data_struct.tissues = [repelem({'AD'},1,length(patient_idxs)),repelem({'NCI'},1,length(NCI))];
    data_struct.genes = geneCounts.gene;
    data_struct.levels  = [normVals(:,patient_idxs), normVals(:,NCI)];
    try
        [clusterModel, ~, addedRxnsForTasks, deletedRxnsInINIT]= ftINIT(prepDataHumanGEM,'AD', [], [], data_struct, {}, getHumanGEMINITSteps('1+1'), true, true);
        clusterModel.id = ['ROSMAP_AD_cluster_' num2str(i)];
        writeYAMLmodel(model,['../../../../models/ROSMAP_AD_cluster_' num2str(i) '.yml'])
        save(['../../../../models/ROSMAP_AD_cluster_' num2str(i) '.mat'],'clusterModel')
    catch
    end
end
