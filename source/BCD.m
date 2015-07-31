%% BCD - block-coordinate descent.
% [Y, iter]  = BCD(S, T, G, Y_prev, M_, lambda, epsilon, k, mu, DEBUG)
% 	S - active set size;
% 	T - number of time samples
% 	G - forward model matrix (in pairs space)
%	Y_prev - solution on previous iteration
% 	M_ - EEG\MEG measurements 
% 	lamdbda - regularization parameter
% 	epsilon - stopping criterion constant
% 	k - number of current IrMxNE iteration (for output)
% 	mu - learning rate
% 	DEBUG - debugging flag 
% 
% 	Y - solution 
% 	iter - total number of iterations made

function [Y, iter]  = BCD(S, T, G, Y_prev, M_, lambda, epsilon, k, mu, DEBUG)
	I = 1000000;			% Number of BCD iterations per one MxNE iteration
	Y_next = Y_prev;
	R = M_ - G * Y_prev;
	fprintf('BCD');
	tic;
	
	for i = 1:I
        range = 1:4;
		for s = 1:S
			Y_next(range,:) = Y_prev(range,:) + mu * G(:,range)' *	R;
			Y_next(range,:) = Y_next(range,:) * max(1 -  mu * lambda / (norm( Y_next(range,:), 'fro' )  ), 0);
			R = R - G(:,range) * ( Y_next(range,:) - Y_prev(range,:) );
            Y_prev = Y_next;
			range = range + 4;
		end
		
		% if DEBUG == true
		% 	fprintf('delta = %g,', norm(Y_next - Y_prev, inf));
		% end

		% if norm(Y_next - Y_prev, inf) < 	% 	fprintf('breaked BCD, delta it = %d, MxNE_it = %d\n', i, k);
		% 	break;
		% end
		% Should compute a dual gap here and check for the convergence
% ------------------------------------------------------------------------------------ %
		
		eta  = dual_gap( [real(M_), imag(M_)], G, [real(Y_next), imag(Y_next)], lambda, S, [real(R), imag(R)] );
% ------------------------------------------------------------------------------------ %
		if i == 1
			fprintf('\nStarting with eta = %f\n', eta);
		end
		if mod(i,5000) == 0 
			fprintf('\n eta = %f\n', eta);
		end
		if mod(i, 10000) == 0
			fprintf('.');			 	 
		end
		if eta < epsilon
			fprintf('\nbreaked BCD, dual it = %d, MxNE_it = %d', i, k);
			break;
		elseif eta > 1e8 && mod(i,1000) == 0
			fprintf('BCD ERROR: diverged!!\n');
			% mu = mu /10;
			% break;
		end
	end
	fprintf('\n');
	fprintf('Done. ');
	toc;
	fprintf('\n');
	fprintf('eta_bcd = %f, ', eta );
	Y = Y_next;
	iter = i;