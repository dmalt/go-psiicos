%% sqToLin: Takes indices i and j as arguments and returnes linear index s
function [s] = sqToLin(i, j, N)
	s = i + (j - 1) * N;
