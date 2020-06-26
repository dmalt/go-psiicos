function [Y, iter]  = BCD(G, Y_prev, M_, lambda, epsilon, mu)
% BCD: block-coordinate descent.
% FORMAT:
%   [Y, iter]  = BCD(G, Y_prev, M_, lambda, epsilon, k, mu)
% INPUTS:
%   G              - {Nsensors^2 x Nsources^2} cross-spectrum generating matrix (squared dimension)
%   Y_prev         - {N_active_pairs x } solution on previous iteration
%   M_             - {Nsensors^2 x Time} EEG\MEG measurements 
%   lamdbda        - regularization parameter
%   epsilon        - stopping criterion constant
%   k              - number of current IrMxNE iteration (for output)
%   mu             - learning rate
% OUTPUTS:
%   Y              - {N_active_pairs x Time} solution matrix 
%   iter           - total number of iterations
% ____________________________________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    tic;
    [Y_,it] = BCD_fast(G, conj(Y_prev'), M_, lambda, epsilon, mu);
    fprintf('\n');
    fprintf('Done. ');
    toc;
    fprintf('\n');
    % fprintf('eta_bcd = %f,\n ', eta );
    Y = conj(Y_');
    iter = it;
end
