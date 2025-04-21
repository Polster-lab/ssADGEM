%load('../results/model_generation.mat')
annotation  = readtable('../data/ROSMAP_annotation_processed.txt','Delimiter','\t');
cluster1 = readtable(['../data/samples_GRN_cluster_1.txt'],'Delimiter',',','NumHeaderLines',0);
cluster2 = readtable(['../data/samples_GRN_cluster_2.txt'],'Delimiter',',','NumHeaderLines',0);
cluster3 = readtable(['../data/samples_GRN_cluster_3.txt'],'Delimiter',',','NumHeaderLines',0);
clusters = [{cluster1} {cluster2} {cluster3}];

cluster1 = strrep(cluster1.Properties.VariableNames,'x','');
cluster2 = strrep(cluster2.Properties.VariableNames,'x','');
cluster3 = strrep(cluster3.Properties.VariableNames,'x','');
annotation.cluster = zeros(height(annotation),1);
[ia,ib] = ismember(cluster1,annotation.patient);
annotation.cluster(ib(ib>0)) = 1;
[ia,ib] = ismember(cluster2,annotation.patient);
annotation.cluster(ib(ib>0)) = 2;
[ia,ib] = ismember(cluster3,annotation.patient);
annotation.cluster(ib(ib>0)) = 3;
%get unique added rxns
allRxns = ihuman.rxns;
tInitAddRxns = getModRxnStats(ihuman,allRxns,annotation,addedRxns_SS,'../results/tINITaddRxns.txt');
tInitRmvRxns = getModRxnStats(ihuman,allRxns,annotation,deletedRxns_SS,'../results/tINITrmvRxns.txt');
C = 94;
D = 70;
fdr_threshold = 0.05;
added = tInitAddRxns(:,[5,9]);
added = added{:,:};
rmved = tInitRmvRxns(:,[5,9]);
rmved = rmved{:,:};

[p_values, bh_corrected_p_added, chi2_stats, is_significant_added] = chi2_feature_enrichment(added, C, D, fdr_threshold);
[p_values, bh_corrected_p_rmved, chi2_stats, is_significant_rmved] = chi2_feature_enrichment(rmved, C, D, fdr_threshold);


function modifiedRxns = getModRxnStats(ihuman,allRxns,annotation,modRxns_SS,filename)
rxnsMat = zeros(length(allRxns),length(modRxns_SS));
for i=1:length(modRxns_SS)
    rxns = modRxns_SS{i};
    if ~isempty(rxns)
        [~,ib] = ismember(rxns,allRxns);
        rxnsMat(ib,i) = 1;
    end
end
sumas = sum(rxnsMat,2);
toKeep = find(sumas>0);
rxnsMat = rxnsMat(toKeep,:);
%all modified Rxns
formulas = constructEquations(ihuman,ihuman.rxns(toKeep));
modifiedRxns = table(ihuman.rxns(toKeep),ihuman.rxnNames(toKeep),formulas,sumas(toKeep));
%now check which models have them and which not (AD vs. NCI)
ADpos = find(strcmpi(annotation.AD,'AD'));
NCIpos = find(strcmpi(annotation.AD,'No_AD'));
ADpresence = zeros(height(modifiedRxns),1);
NCIpresence = zeros(height(modifiedRxns),1);
c1 = zeros(height(modifiedRxns),1);
c2 = zeros(height(modifiedRxns),1);
c3 = zeros(height(modifiedRxns),1);
for i=1:height(modifiedRxns)
    idxs = find(rxnsMat(i,:));
    [ia,ib_AD] = ismember(idxs,ADpos);
    ADpresence(i) = sum(ia);
    c1(i) = sum(annotation.cluster(ia)==1);
    c2(i) = sum(annotation.cluster(ia)==2);
    c3(i) = sum(annotation.cluster(ia)==3);
    [ia,ib_NCI] = ismember(idxs,NCIpos);
    NCIpresence(i) = sum(ia);   
end
modifiedRxns.ADpresence = ADpresence;
modifiedRxns.c1 = c1/38;
modifiedRxns.c2 = c2/33;
modifiedRxns.c3 = c3/20;
modifiedRxns.NCIpresence = NCIpresence;
modifiedRxns.ADperc = modifiedRxns.ADpresence/numel(ADpos);
modifiedRxns.NCIperc = modifiedRxns.NCIpresence/numel(NCIpos);
modifiedRxns.presenceRatio = modifiedRxns.ADperc./modifiedRxns.NCIperc;
modifiedRxns = sortrows(modifiedRxns,{'ADperc' 'presenceRatio'},{'descend' 'descend'});
writetable(modifiedRxns,filename,'Delimiter','\t');
end