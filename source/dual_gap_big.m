function eta = dual_gap(M, G, X, lambda, S, R, k, w)
	% B = sqrt(max(diag(G' * R * R' * G)))/ lambda;
	B(S) = 0;
	sigma = 0;
	range = 1:4;
	if k ==1
		B = load('../aux/V.txt');
	elseif k ~= 1
		idx = support(w);
		[dummy, Sr] = size(idx);
		B(S) = 0; % Memory allocation. 
		for sr = 1:Sr
			B(sr)  = norm(G_pair(idx(sr),G, w(idx(sr)))' * R, 'fro');
		end	
	end
	for s=1:S
		sigma = sigma + norm(X(range,:), 'fro');
		range = range + 4;
	end
	S
	[C, I] = max(B)
	C = C / lambda
	R_ = R / max(C, 1);
	At = 0.5 * (norm(R, 'fro') ^ 2 )
	Bt = lambda * sigma
	Ct = 0.5 * ( norm(R_,'fro') ^ 2 )
	Dt = sum(sum(R_ .* M ))
	eta = At + Bt + Ct - Dt;
	