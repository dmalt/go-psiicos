%% dual_gap: returns primal-dual gap
function eta = dual_gap(M, G, X, lambda, S, R)
    % B = sqrt(max(diag(G' * R * R' * G)))/ lambda;
    B = zeros(S,1);
    sigma = 0;
    range = 1:4;
    for s=1:S
        B(s) = norm(G(:,range)' * R, 'fro');
        sigma = sigma + norm(X(range,:), 'fro');
        range = range + 4;
    end
    C = max(B);
    C = C / lambda;
    R_ = R / max(C, 1);
    eta = 0.5 * (norm(R, 'fro') ^ 2 ) + lambda * sigma + 0.5 * ( norm(R_,'fro') ^ 2 ) - sum(sum(R_ .* M ));
