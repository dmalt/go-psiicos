function supp = support(X)
%% support: returns nonzero raws of a matrix
% FORMAT:
%    supp = support(X)
% INPUTS:
%   X            - {N x M} matrix
% OUTPUTS:
%   supp         - {Nnonzero x M} matrix of nonzero raws of X
% ___________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    [S,T] = size(X);
    supp = [];
    for s = 1:S
        if norm(X(s,:) , inf) > 0
            supp = [supp, s];
        end
    end
end
