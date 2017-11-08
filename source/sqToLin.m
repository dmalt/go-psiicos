function [s] = sqToLin(i, j, N)
% -----------------------------------------------------------------------
% sqToLin: Takes indices i and j as arguments and returnes linear index s
% -----------------------------------------------------------------------
% description
% FORMAT:
%   format 
% INPUTS:
%   inputs        -
% OUTPUTS:
%   outputs
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    s = i + (j - 1) * N;
end
