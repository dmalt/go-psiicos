%% linToSq.m: Takes linear index s as an argument. Returns indices i and j				
function [i, j] = linToSq(s, N)
	i = mod(s, N);
	if i == 0 
		i = N;
	end
	j = (s - i) / N + 1;
