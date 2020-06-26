function A = ActiveSet(G_small, R_, lambda, w, k, V)
%% ActiveSet: compute an array of indices of sources which violate KKT conditions the most
% FORMAT:
%   A = ActiveSet(G_small, R_, lambda, w, k, V)
% INPUTS:
%   G_small         - {Nsensors x Nsources} forward model matrix
%   R_              - {}
%   lambda          - regularization parameter. 
%   w               - 
%   V               - active set chunck
% OUTPUTS:
%   A               - indeces of V first sources which violate KKT conditions the most
% _________________________________________________________________________________________
% Dmitrii Altukhov, dm.altukhov@me.com

    fprintf('Calculating active set...\n');
    fprintf('Start (Active set)\n');
    startBCD = tic;
    R = [real(R_), imag(R_)];
    if k == 1       
        idx = @(x) x;
        S = length(w);
        save ../aux/G_small.txt G_small -ASCII
        save ../aux/R.txt R -ASCII
        unix('./ActSetViol');
        src_violations = load('../aux/V.txt');
        fprintf('Stop (Active set). ');
        src_violations = src_violations' - lambda;
    elseif k ~= 1
        idx = support(w);
        [dummy, S] = size(idx);
        src_violations(S) = 0; % Memory allocation. 
        for s = 1:S
            src_violations(s)  = norm(G_pair(idx(s), G_small, w(idx(s)))' * R, 'fro') - lambda;
        end 
    end
    [Dummy, A] =  sort(src_violations, 'descend');

    if V < S
         A = A(1:V);
    end

    
    % src_violations(1:25)
    A = A(src_violations(A) > 0);
    % src_violations(A)
    A = idx(A);
    toc(startBCD);
    fprintf('Done (Active set)\n');
    A = sort(A);
end
