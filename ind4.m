%% ind4: transforms indices from pairs of sites to pairs of sources		
function [idx4] = ind4(idx)
	S = length(idx);
	idx4 = [];
	for s=1:S
		idx4 = [idx4,idx(s)*4-3:idx(s)*4];
	end
