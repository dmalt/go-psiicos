function [idx4] = ind4(idx)
% ----------------------------------------------------------------
% ind4: transforms indices from pairs of sites to pairs of sources      
% ----------------------------------------------------------------
% FORMAT:
%   idx4 = ind4(idx) 
% INPUTS:
%   idx        - index of pair of sites
% OUTPUTS:
%   idx4       - indeces of pairs of interacting
%                dipoles for site idx
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    S = length(idx);
    idx4 = [];
    for s=1:S
        idx4 = [idx4,idx(s)*4-3:idx(s)*4];
    end
end
