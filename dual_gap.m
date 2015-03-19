%% dual_gap: returns primal-dual gap
function eta = dual_gap(M, G, X, lambda, S, R)
	R_ = R / max(sqrt(max(diag(G' * R * R' * G)))/ lambda, 1);
	eta = 0.5 * (norm(M - G * X, 'fro') ^ 2 )+ lambda * sum(sqrt(diag(X * X'))) + 0.5 * ( norm(R_,'fro') ^ 2 ) - sum(sum(R_ .* M ));