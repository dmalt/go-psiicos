%% G_pair: returns topograhy matrix of interaction of i and j sites
function G = G_pair(s, G_small, w)
	G = [];
	[Nch, Nsrc] = size(G_small);
	Nsites = Nsrc / 2;
	[i, j] = linToSq(s, Nsites); % i-th and j-th dipoles are interacting
	s11 = (s - 1) * 4 + 1;
	s12 = (s - 1) * 4 + 2;
	s21 = (s - 1) * 4 + 3;
	s22 = (s - 1) * 4 + 4;
	G = [columnG_fast(s11, G_small, w), columnG_fast(s12, G_small, w), columnG_fast(s21, G_small, w), columnG_fast(s22, G_small, w)];
