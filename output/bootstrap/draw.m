% ISD = load('InputData4Simulations.mat');
Ctx = load('/home/meg/fif_matlab/Brainstorm_db/PSIICOS/anat/0003_pran/tess_cortex_concat_2000V.mat');
FolderName = '/home/meg/fif_matlab/Brainstorm_db/PSIICOS/data/';
ChUsed = 1:306; ChUsed(3:3:end) = [];
TimeRange = [0, 0.700];
Subject = '0003_pran/brainstormsubject.mat';
%Subject = '0019_shev/brainstormsubject.mat';

Conditions = {'1','4'}; % '2','4'};
Band = [16 25];
%Band = [70 85];
% Band = [8 12];
Fsamp = 500;
[b,a] = butter(5,Band/(Fsamp/2));
Protocol = bst_get('ProtocolStudies','PSIICOS');
clear ConData;
fprintf('Loading real data from BST database.. \n');
%% Load data and compute cross-spectral matrix 
ConditionsFound = 0;
clear ConData;
for c = 1:length(Conditions)
    for s = 1:length(Protocol.Study)
        if(strcmp(Protocol.Study(s).Name,Conditions{c}) & strcmp(Protocol.Study(s).BrainStormSubject,Subject))
            fprintf('Found study condition %s %s\n ', Conditions{c},Protocol.Study(s).BrainStormSubject); 
            for hm = 1:length(Protocol.Study(s).HeadModel)
                if(strcmp(Protocol.Study(s).HeadModel(hm).Comment,'Overlapping spheres_HR'))
                    ConData{c}.HM_HR = load([FolderName Protocol.Study(s).HeadModel(hm).FileName]);
                else
                    ConData{c}.HM_LR = load([FolderName Protocol.Study(s).HeadModel(hm).FileName]);
                end
            end;
        end;
    end;
end;

for c = 1:length(Conditions)
    for s = 1:length(Protocol.Study)
        if(strcmp(Protocol.Study(s).Name,Conditions{c}) & strcmp(Protocol.Study(s).BrainStormSubject,Subject))
            fprintf('Found study condition %s %s\n ', Conditions{c},Protocol.Study(s).BrainStormSubject); 
            ConData{c}.NumTrials = length(Protocol.Study(s).Data);
            fprintf('Loading Trials (Max %d) : ', ConData{c}.NumTrials); 
            for t = 1:ConData{c}.NumTrials
                aux = load([FolderName Protocol.Study(s).Data(t).FileName]);
                if(t==1)
                    [ans, ind0] =min(abs(aux.Time-TimeRange(1)));
                    [ans, ind1] =min(abs(aux.Time-TimeRange(2)));
                    T = ind1-ind0+1; 
                    ConData{c}.Trials = zeros(Nch,T,ConData{c}.NumTrials);
                    ConData{c}.Fsamp = 1./(aux.Time(2)-aux.Time(1));
                end;
                tmp = filtfilt(b,a,(UP*aux.F(ChUsed,:))')';
                ConData{c}.Trials(:,:,t) = tmp(:,ind0:ind1);
                if t>1
                    for tt=0:log10(t-1)
                        fprintf('\b'); % delete previous counter display
                    end
                end
                fprintf('%d', t);
            end; % trials t
            fprintf(' -> Done\n');
        end;
    end;
    
    if(length(ConData)>=c)
        P = sum(sum(abs(ConData{c}.Trials),1),2);
        Pm = median(squeeze(P));
        ind = find(P>2*Pm | P<0.5*Pm);
        REJ{c} = ind;
        ConData{c}.Trials(:,:,ind) = [];
        fprintf('Computing cross-spectral matrix ....' ); 
        ConData{c}.CrossSpec = CrossSpectralMatrix(ConData{c}.Trials,Band,500);
        ConData{c}.CrossSpecTime = CrossSpectralTimeseries( ConData{c}.Trials);
        %ConData{c}.CrossSpec = reshape(mean(ConData{c}.CrossSpecTime,2),Nch,Nch);
        fprintf('-> Done\n' ); 
    end;
end;

figure;
hctx  = trisurf(Ctx.Faces,Ctx.Vertices(:,1),Ctx.Vertices(:,2),Ctx.Vertices(:,3),'FaceColor',[0.1,0.51,1], 'EdgeColor','none','FaceAlpha', 0.1);
hold on;

camlight left; lighting phong
camlight right; 
hold on;
 % plot3(XYZGen(:,1), XYZGen(:,2),XYZGen(:,3),'r');
cols = ['g','r','m','y','k','b'];

R = ConData{1}.HM_LR.GridLoc;%GLowRes.GridLoc;
Dmax = 0.02;
% NetworkPairIndex{2} = [1,2,3];
NPI = [1, 2, 3];  %NetworkPairIndex{2};
% create binary arrays indicators for each network from NPI 
Nsites = size(R, 1);
A = [];
out = dir('./Output_*');
numfiles = length(out);
mydata = cell(1, numfiles);
for k=1:numfiles
	mydata{k} = load(out(k).name);
	A = [A,mydata{k}.A];
end
% B = setdiff(A, unique(A));
occ = zeros(size(A));
for i = 1:length(A)
	occ(i) = sum(A == A(i));
end
B = A(occ>0);
Result = [(mod(B,Nsites))', ((B - mod(B,Nsites)) / Nsites+ 1)'];
[Npairs, dummy] = size(Result);
adjMat = zeros(Npairs,Npairs);
Dpair = 0.013;
for p1 = 1:Npairs
    for p2 = p1:Npairs
        if norm(R(Result(p1,1),:) - R(Result(p2,1),:)) < Dpair && norm(R(Result(p1,2),:) - R(Result(p2,2),:)) < Dpair
            adjMat(p1,p2) = 1;
            adjMat(p2,p1) = 1;
        end
    end
end
Result ;
N = size(Result,1);
i = 1;
while(~isempty(Result))
    if length(nonzeros(bfs(adjMat,i) > 0)) > 10
        % Result(i,:)
        clust = Result(bfs(adjMat,i) > -1,:);
        drawset(clust, R);
        
        Result = Result(bfs(adjMat,i) == -1,:);
        adjMat = adjMat(bfs(adjMat,i) == -1, bfs(adjMat,i) == -1);
    	% L1 = line( R(Result(i,:), 1), R(Result(i,:),2), R(Result(i,:),3) );
    	% plot3( R(Result(i,1), 1), R(Result(i,1),2), R(Result(i,1),3), 'c.', 'Markersize', 40 );
     %    plot3( R(Result(i,2), 1), R(Result(i,2),2), R(Result(i,2),3), 'c.', 'Markersize', 40 );
    	% set(L1, 'Color', 'c', 'linewidth',2);
        % i = 1;
    else
        Result = Result(2:end,:);
        adjMat = adjMat(2:end,2:end);
        % i = i + 1;
    end
end
hold on;
collim = [0,15e-3];

CT2 = ConData{1}.CrossSpecTime;
ncomp = 5;
X_res = cell(1,numfiles)
for k=1:numfiles
    CT = mydata{k}.BootsCT;
    M  = ProjOut(CT, CT2, G2dLRU) ;
    M_abs = M / norm(M);
    
    [Mu Ms Mv] = svd(M_abs);
    X_ = mydata{k}.X;
    X_res{k} = X_ * Mv(:,1:ncomp)';
end


