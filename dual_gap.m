%% dual_gap: returns primal-dual gap
function [eta, Fp, Fd] = dual_gap(M, G, X, lambda, S, R)
	Sigma_groups = 0;
	tmp = zeros(S,1);

	for s = 1:S
		Sigma_groups = Sigma_groups + norm(X(s,:), 2);
		tmp(s) = norm(G(:,s)' * R, 'fro');
	end

	Fp = 0.5 * (norm(M - G * X, 'fro') ^ 2 )+ lambda * Sigma_groups;
	R_ = R / max(max(tmp) / lambda, 1);
    R_;
	% R_ = R / max(tmp) / lambda;
	Fd =  -0.5 * norm(R_,'fro')  + trace(R_' * M);
	eta = Fp - Fd;