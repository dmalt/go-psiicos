%% CalcG: returns subset of columns of matrix G, which is really helpful, because G is super big. idx - vector of column numbers
function [subG] = CalcG(idx, G_small, w)
    subG = [];
    for i = idx
        subG = [subG, G_pair(i, G_small, w(i))];
    end
