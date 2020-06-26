% ------------------------------------------ %
% ------------------------------------------ %
% DATE: 2020-06-14
% AUTHOR: dmalt
% ------------------------------------------ %

% ---------------- simulate data --------------- %
subjID = 'test';
% PhaseLag = pi / 4 + pi / 2;
PhaseLag = pi / 2 - pi / 20;
% PhaseLag = pi / 4;
% PhaseLag = pi / 4;
% GainSVDTh = 0.001;
GainSVDTh = 0.01;
% GainSVDTh = 0;
n_tr = 100;
snr_induced = 0.5;
snr_evoked = 0;
is_induced = false;
jitter = 0.1;

[HM, Trials, Ctx] = ups.sim.SimulateData( ...
    PhaseLag, n_tr, GainSVDTh, snr_induced, snr_evoked, jitter);
CT = ups.conn.CrossSpectralTimeseries(Trials, true);
pwr_rnk = 500;
CT_proj = ps.ProjectFromSlComplete(CT, HM.gain, pwr_rnk);

lambda = 0.1;
ActSetChunk = 100;
t1 = tic;
[X,Aidx] = go_psiicos(CT, [], HM.gain, lambda, ActSetChunk, 'sim_3ntw');
t2 = toc(t1);
fprintf('Total time: %f\n', t2)

clear con_idx
k = 1;
for i_pair = 1:2:length(Aidx)
    [i,j] = linToSq(Aidx(i_pair), length(HM.GridLoc));
    con_idx{k} = [i, j];
    k = k + 1;
end

con = ups.Bundles(con_idx, HM, Ctx);

con.PlotViews(0.05, 2, 0.004);
set(gcf, 'Color', 'w')

return 

colors = ups.plt.GetColors(3)
c(1,:) = colors(2).RGB
c(2,:) = colors(3).RGB
c(3,:) = colors(4).RGB

xx = reshape(X', 500, 4, []);
xxm =(squeeze(mean(abs(xx), 2)));

figure;
t = linspace(0, 1000, length(xxm(:,1)));
plot(t, xxm(:,1), 'Color', c(1,:), 'LineWidth', 2)
hold on;
plot(t, xxm(:,3), 'Color', c(2,:), 'LineWidth', 2)
hold on;
plot(t, xxm(:,5), 'Color', c(3,:), 'LineWidth', 2)
set(gca, 'ColorOrder', c, 'FontSize', 20)
grid 
set(gcf, 'Color', 'w')
xlabel('Time, ms')

export_fig('timeseries_gopsiicos_pi2pi20', '-jpg')
export_fig('networks_gopsiicos_pi2pi20', '-jpg')
% load handel.mat;
% yRange = [-0.7,0.7];
% soundsc(y,yRange);
