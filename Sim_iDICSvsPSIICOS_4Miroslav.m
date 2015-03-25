clear all;
close all;
%% Parameters block
InducedScale = 0.35;
EvokedScale = 0.0/2.5;
GainSVDTh = 0.001; % 0.05 results into 47 eigensensors and makes it run faster but produces less contrasting subcorr scans
% for a more reliable preformance use 0.01 to get all the sensor on board but be ready to wait;
NetworkPairIndex{1} = [1,2];
NetworkPairIndex{2} = [1,2,3];
%% Load forward model and reduce it  
% load reduced forward model (GLowRes)
load('GLowRes.mat'); 
% get grid node locations
R = GLowRes.GridLoc;
% set to use gradiometers only
ChUsed = 1:306; 
ChUsed(3:3:end) = [];
% calculate tangential plane dipoles
[Nch, Nsites] = size(GLowRes.Gain(ChUsed,1:3:end));
G2d = zeros(Nch,Nsites*2);
G2d0 = zeros(Nch,Nsites*2);
range = 1:2;
for i=1:Nsites
    g = [GLowRes.Gain(ChUsed,1+3*(i-1)) GLowRes.Gain(ChUsed,2+3*(i-1)) GLowRes.Gain(ChUsed,3+3*(i-1))];
    [u sv v] = svd(g);
    gt = g*v(:,1:2);
      G2d(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
      G2d0(:,range) = gt;
    range = range + 2;
end;
% reduce the sensor space
[ug sg vg] = spm_svd(G2d*G2d',GainSVDTh);
UP = ug';
G2dU = UP*G2d;
G2d0U = UP*G2d0;

PHI = [pi/20 pi/2];
it = 1;
for phi = PHI 
    %% Data simulation
    if(it==1)
        % we will use the same phase shifts in the second and subsequent
        % iterations. We will store the random phases in PhaseShifts 
        [ Evoked, Induced, BrainNoise, SensorNoise, G2dHiRes, RHiRes, Fs,Ntr, XYZGen, Ggen, PhaseShifts] = ...
        SimulateDataPhase(NetworkPairIndex{2},phi,true,[]);
    else
        % and use PhaseShits from the first iteration
        [ Evoked, Induced, BrainNoise, SensorNoise, G2dHiRes, RHiRes, Fs,Ntr, XYZGen, Ggen, PhaseShifts] = ...
        SimulateDataPhase(NetworkPairIndex{2},phi,false,PhaseShifts);
    end;        
    % mix noise and data 
    % in order to control SNR we first normalize norm(BrainNoise(:)) = 1 and 
    % norm(Induced(:)) = 1 and then mix the two with the coefficient
    Data0 = BrainNoise + InducedScale*Induced + EvokedScale*Evoked ;
    [bf af] = butter(5,[8,12]/(0.5*500),'bandpass');
    % Filter in the band of interest
    Data = filtfilt(bf,af,Data0')';
    %% Reshape the data in a 3D structure(Space x Time x Epochs)
    [Nch Tcnt] = size(Data);
    T = fix(Tcnt/Ntr);
    Nch = size(UP,1);
    %reshape Data and store in a 3D array X
    X1 = zeros(Nch,T,Ntr);
    range = 1:T;
    for i=1:Ntr
        X1(:,:,i) = UP*Data(:,range);
        range = range+T;
    end;
    %% Calculate band cross-spectral matrix 
    CrossSpecTime{it} = CrossSpectralTimeseries(X1);
    C{it} = reshape(mean(CrossSpecTime{it},2),Nch,Nch); 
    
    % %% Experiment with different methods
    [Qpsiicos{it}, IND{it}, Cp{it} ] = PSIICOS( C{it}, G2dU);
    it = it+1;
end;
return;
