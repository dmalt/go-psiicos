function  bootstGoPs(Subj)

%% bootstGoPs: Launches bootstrap for a subject Subj
% function bootstGoPs(Subj)
path(path, '../input/');
[sSubjData, G2dLRU] = PrepRealData(Subj);
BootsTrials = zeros(size(sSubjData{2}.Trials));

for it=1:100
    redNumTr = size(sSubjData{2}.Trials,3);
    resample = randi(redNumTr, 1, redNumTr);
    for tr = 1:redNumTr
        BootsTrials(:,:,tr) =  sSubjData{2}.Trials(:,:,resample(tr));
    end
    BootsCT = CrossSpectralTimeseries( BootsTrials);
    cd ../source/;
    [X,A] = go_psiicos(0.08, 500, G2dLRU, BootsCT, sSubjData{1}.CrossSpecTime, Subj);
    cd ../input/;
    it_str = num2str(it);
    dirName = strcat('bootstr_',Subj,'/');
    mkdir('../output', dirName);
    outDir = strcat('../output/', dirName);
    save ( strcat( strcat(outDir,'Output_', it_str), '.mat'), 'A', 'X', 'resample', 'BootsCT');
end
