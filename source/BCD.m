%% BCD - block-coordinate descent.
% [Y, iter]  = BCD(S, T, G, Y_prev, M_, lambda, epsilon, k, mu)
% 	S - active set size;
% 	T - number of time samples
% 	G - forward model matrix (in pairs space)
%	Y_prev - solution on previous iteration
% 	M_ - EEG\MEG measurements 
% 	lamdbda - regularization parameter
% 	epsilon - stopping criterion constant
% 	k - number of current IrMxNE iteration (for output)
% 	mu - learning rate
% 
% 	Y - solution 
% 	iter - total number of iterations made

function [Y, iter]  = BCD(G, Y_prev, M_, lambda, epsilon, mu)
	tic;

	[Y_,it] = BCD_fast(G, conj(Y_prev'), M_, lambda, epsilon, mu);
	fprintf('\n');
	fprintf('Done. ');
	toc;
	fprintf('\n');
	% fprintf('eta_bcd = %f\n, ', eta );
	Y = conj(Y_');
	iter = it;

