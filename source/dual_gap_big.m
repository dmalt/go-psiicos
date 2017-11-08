function eta = dual_gap(M, G, X, lambda, S, R)
%% Calculate duality gap for M, G, X, lambda, S, R
% 
% FORMAT:
%   eta = dual_gap(M, G, X, lambda, S, R)
% INPUTS:
%   M          - {Nsources^2 x Time} measured cross-spectrum
%   G          - {Nsensors^2 x Nsources^2} cross-spectrum generating matrix for pairs of connections
%   X          - {Nsources^2 x Time} current solution timecources
%   lambda
%   S          - number of pairs of sources; equal or less then the total number of 
%                pairs of sources on cortex
%   R          - {Nsensors^2 x Time} - residual error; R = M - G * X 
% OUTPUTS:
%   eta        - dual gap
% ________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
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
    eta = 0.5 * (norm(R, 'fro') ^ 2 ) + lambda * sigma + 0.5 * ( norm(R_,'fro') ^ 2 ) - sum(sum(R_ .* M ));
end    
