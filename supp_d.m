%% supp_d: calculates support set for given matrix X each 4 rows of wich represent single dipole pair
function [supp] = supp_d(X)
	[Nsrc_pairs, T] = size(X);
	Nsite_pairs = Nsrc_pairs / 4;
	range = 1:4;
	normX(Nsite_pairs) = 0;
	for s = 1:Nsite_pairs
		normX(s) = norm(X(range,:), 'fro');
		range =range + 4;
	end
	supp = support(normX');
