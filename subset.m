%% subset: returns an array of indices of sources which violate KKT conditions the most
function A = subset(G, R, lambda, S)
	V = 20; % Number of sources in a subset
	src_violations = zeros(1,S);
	for s = 1:S
		src_violations(s)  = norm(G(:,s)' * R, 2) - lambda;
	end	
	[Dummy, A] =  sort(src_violations, 'descend');
	A = A(1:V);
