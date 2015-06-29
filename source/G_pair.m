%% G_pair: returns topograhy matrix of interaction of i and j sites
% s is number of a pair of sites
% G_small is generating model matrix in nonsquare space
% w is scalar weight of pair s
function G = G_pair(s, G_small, w)
	G = [];
	[Nch, Nsrc] = size(G_small);
	Nsites = Nsrc / 2;
	[i, j] = linToSq(s, Nsites); % i-th and j-th dipoles are interacting
	s11 = sqToLin(2 * i - 1, 2 * j - 1, Nsrc);
	s12 = sqToLin(2 * i, 2 * j - 1, Nsrc);
	s21 = sqToLin(2 * i - 1, 2 * j, Nsrc);
	s22 = sqToLin(2 * i, 2 * j, Nsrc);
	G = [columnG_fast(s11, G_small, w), columnG_fast(s12, G_small, w), columnG_fast(s21, G_small, w), columnG_fast(s22, G_small, w)];
