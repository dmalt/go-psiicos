%% support: returns nonzero strings of a matrix
function supp = support(X)
	[S,T] = size(X);
	supp = [];
	for s = 1:S
		if norm(X(s,:) , inf) > 0
			supp = [supp, s];
		end
	end
