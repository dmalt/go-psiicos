%% dual_gap: returns primal-dual gap
function [eta, Fp, Fd] = dual_gap(M, G, X, lambda, S, R)
	Sigma_groups = 0;
	tmp(S) = 0;

	for s = 1:S
		Sigma_groups = Sigma_groups + norm(X(s,:), 2);
		tmp(s) = norm(G(:,s)' * R);
	end

	Fp = 0.5 * (norm(M - G * X, 'fro') ^ 2 )+ lambda * Sigma_groups;
	R_ = R / max(max(tmp) / lambda, 1);
	% R_ = R / max(tmp) / lambda;
	Fd1 = 0.5 * ( norm(R_,'fro') ^ 2 );
	Fd2 = sum(sum(R_ .* M ));
	Fd =  Fd2 - Fd1;
	eta = Fp - Fd;