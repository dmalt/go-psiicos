function G = G_pair(s, G_small, w)
%% G_pair: returns topograhy matrix of interaction of i and j sites
%          multiplied by weight w
% FORMAT:
%   G = G_pair(s, G_small, w) 
% INPUTS:
%   s               - index of a pair of interacting sites (s = i + Nsources * j)
%   G_small         - {Nsources x Nsensors} forward model matrix
%   w               - weight for s dipole pair
% OUTPUTS:
%   G               - {Nsensors ^ 2 x 4}
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    G = [];
    [Nch, Nsrc] = size(G_small);
    Nsites = Nsrc / 2;
    [i, j] = linToSq(s, Nsites); % i-th and j-th dipoles are interacting
    s11 = (s - 1) * 4 + 1;
    s12 = (s - 1) * 4 + 2;
    s21 = (s - 1) * 4 + 3;
    s22 = (s - 1) * 4 + 4;
    G = [columnG_fast(s11, G_small, w), columnG_fast(s12, G_small, w), columnG_fast(s21, G_small, w), columnG_fast(s22, G_small, w)];
end
