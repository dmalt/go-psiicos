function [i, j] = linToSq(s, N)
% ----------------------------------------------------------------------
% linToSq.m: takes linear index s and liner dimension N as an arguments,
% returns indices i and j; s = i + N * (j - 1)
% ----------------------------------------------------------------------
% FORMAT:
%   [i, j] = linToSq(s, N) 
% INPUTS:
%   s         - linear index 
%   N         - linear dimension 
% OUTPUTS:
%   i         - first index
%   j         - second index
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    i = mod(s, N);
    if i == 0 
        i = N;
    end
    j = (s - i) / N + 1;
end
