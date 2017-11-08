function [supp] = supp_d(X)
%% supp_d: calculates support set for given matrix X each 4 rows 
%          of which represent single dipole pair; set of 4 raws 
%          will be included in supp <=> all 4 raws for this pair 
%          are nonzero
% FORMAT:
%    supp = supp_d(X)
% INPUTS:
%   X          - {4 * Nsource_pairs x Time} matrix of timecources for interacting, some of which
%                are probably zero; timecources go in sets of 4 rows (timecources) for each pairs.
% OUTPUTS:
%   supp       - {4 * Nnonzero_source_pairs x Time} matrix of nonzero timecources
% __________________________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    [Nsrc_pairs, T] = size(X);
    Nsite_pairs = Nsrc_pairs / 4;
    range = 1:4;
    normX(Nsite_pairs) = 0;
    for s = 1:Nsite_pairs
        normX(s) = norm(X(range,:), 'fro');
        range =range + 4;
    end
    supp = support(normX');
end
