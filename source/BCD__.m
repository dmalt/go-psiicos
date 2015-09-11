%% BCD - block-coordinate descent.
% [Y, iter]  = BCD(G, Y_prev, M_, lambda, epsilon, k, mu)
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

function [Y, iter]  = BCD__(G, Y_prev, M_, lambda, epsilon, mu)
	I = 2000000;			% Number of BCD iterations per one MxNE iteration
	Y_next = Y_prev;
	R = M_ - G * Y_prev;
	T = size(M_,2);
	S = size(G,2) / 4;
	eta = 2;
	fprintf('BCD');
	tic;
	% imag(R)
	% imag(M_)
	% ouy =imag(Y_prev);
	for i = 1:I
		s =  randi(S,1,1); %;mod(i,S)%
        range = 4 * s - 3:4 * s;
		Y_next(range,:) = Y_next(range,:) + mu * G(:,range)' *	R;
		Y_next(range,:) = Y_next(range,:) * max(1 -  mu * lambda / (norm( Y_next(range,:), 'fro' )  ), 0);
		R = R - G(:,range) * ( Y_next(range,:) - Y_prev(range,:) );
		Y_prev = Y_next;
% ------------------------------------------------------------------------------------ %
		% eta  = dual_gap( [real(M_), imag(M_)], G, [real(Y_next), imag(Y_next)], lambda, S, [real(R), imag(R)] );
% ------------------------------------------------------------------------------------ %
		if mod(i, 20 * S) == 0
			eta  = dual_gap( [real(M_), imag(M_)], G, [real(Y_next), imag(Y_next)], lambda, S, [real(R), imag(R)] );
			% eta  = dual_gap_fast( [real(M_), imag(M_)], G, [real(Y_next), imag(Y_next)]', lambda, S, [real(R), imag(R)] );
		end
		if mod(i, 10000) == 0
			fprintf('.');			 	 
		end
		if eta < epsilon
		 	fprintf('\nbreaked BCD, dual it = %d', i);
		 	break;
		elseif eta > 1e8 && mod(i,1000) == 0
			fprintf('BCD ERROR: diverged!!\n');
			mu = mu /10;
			% break;
		end
	end
	fprintf('\n');
	fprintf('Done. ');
	toc;
	fprintf('\n');
	fprintf('eta_bcd = %6.6f\n, ', eta );
	Y = Y_next;
	iter = i;
