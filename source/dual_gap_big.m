function eta = dual_gap(M, G, X, lambda, S, R)
	% B = sqrt(max(diag(G' * R * R' * G)))/ lambda;
	B(S) = 0;
	sigma = 0;
	range = 1:4;
	B = load('../aux/V.txt');
	for s=1:S
		sigma = sigma + norm(X(range,:), 'fro');
		range = range + 4;
	end
	C = max(B);
	C = C / lambda;
	R_ = R / max(C, 1);
	fprintf('norm R^2 = %f\n',0.5 * (norm(R, 'fro') ^ 2 ) );
	eta = 0.5 * (norm(R, 'fro') ^ 2 ) + lambda * sigma + 0.5 * ( norm(R_,'fro') ^ 2 ) - sum(sum(R_ .* M ));
	