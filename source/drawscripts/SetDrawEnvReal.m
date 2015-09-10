
path(path,'../../input/');
% Subj = '0019_shev';
Subj = '0003_pran';
[ConData, G2dLRU] = PrepRealData(Subj);
brainstormDir = '~/fif_matlab/Brainstorm_db/';
% anatPath = strcat(brainstormDir,'PSIICOS/anat/', Subj,'/tess_cortex_pial_low_2000V.mat');
anatPath = strcat(brainstormDir,'PSIICOS/anat/', Subj,'/tess_cortex_concat_2000V.mat');
Ctx = load(anatPath);
gopsiicosFolder = '~/Dropbox/Documents/Education/MEG/Osadchii/gopsiicos/gopsiicos_source'


R = ConData{1}.HM_LR.GridLoc;%GLowRes.GridLoc;
% NetworkPairIndex{2} = [1,2,3];
NPI = [1, 2, 3];  %NetworkPairIndex{2};
% create binary arrays indicators for each network from NPI 


% Load bootstrap data
Nsites = size(R, 1);
A = [];
subjBootstrDir = strcat(gopsiicosFolder, '/output/','bootstr_',Subj,'/');
out = dir(strcat(subjBootstrDir, 'Output_*'));
numfiles = length(out);
mydata = cell(1, numfiles);
for k=1:numfiles
	mydata{k} = load(strcat(subjBootstrDir, out(k).name));
    if length(mydata{k}.A) < 100
    	A = [A,mydata{k}.A];
    end
end
% B = setdiff(A, unique(A));


collim = [0,15e-3];
% return;
CT2 = ConData{1}.CrossSpecTime;
ncomp = 5;
X_res = cell(1,numfiles);
M_res = cell(1,numfiles);
for k=1:numfiles
    CT = mydata{k}.BootsCT;
    % M  = ProjOut(CT, CT2, G2dLRU) ;
    % M_abs = M / norm(M);
    
    % [Mu Ms Mv] = svd(M_abs);
    X_ = mydata{k}.X;
    X_res{k} = X_ ; %Mv(:,1:ncomp)';
    M_res{k} = CT;
end

