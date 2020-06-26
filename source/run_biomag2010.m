home = getenv('HOME');
bst_path = [home,'/Documents/MATLAB/bst/brainstorm_db'];
data_path = [home, '/Data/mentrot/MentalRotationDeLange/preprocessed'];
subj_ID = 'biomag2010';
protocol_path = [bst_path,'/mentrot'];
condition = 'raw';

% head model
isLR = true;
GainSVDTh = 0.0001; % results in 45 components
ch_type = 'MEG';


time_range = [0, 1];

alpha_band = [8,12];
beta_band = [16,24];

freq_band = beta_band;
path(path, '/home/dmalt/Code/matlab/PSIICOS/scripts/biomag2010');
% load head model
HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh, ch_type);
ltr     = memoize(@load_trials);
calc_CT = memoize(@ups.conn.CrossSpectralTimeseries);
proj    = memoize(@ps.ProjectFromSlComplete);
msvd    = memoize(@svd);

tr_filt= ltr(data_path, time_range, freq_band, HM);
CT = calc_CT(tr_filt, true);

pwr_rnk = 150;
CT_proj = proj(CT, HM.gain, pwr_rnk);


lambda = 0.05;
ActSetChunk = 40;
[X,Aidx] = go_psiicos(CT, [], HM.gain, lambda, ActSetChunk, 'biomag');


for i_pair = 1:length(Aidx)
    [i,j] = linToSq(Aidx(i_pair), length(HM.GridLoc));
    con_idx{i_pair} = [i,j];
end

[Ctx, CtxHR, CtxHHR] = ups.bst.GetCtx(subj_ID, protocol_path);
con = ups.Bundles(con_idx, HM, CtxHHR)
