function [Y, iter]  = BCD(S, T, G, Y_prev, M, lambda, epsilon, k, mu, tau, DEBUG)
	I = 20000;			% Number of BCD iterations per one MxNE iteration
	Y_next = Y_prev;
	R = M - G * Y_prev;

	for i = 1:I
		for s = 1:S
			Y_next(s,:) = Y_prev(s,:) + mu .* G(:,s)' *	R;
			Y_next(s,:) = Y_next(s,:) * max(1 -  mu * lambda / (norm( Y_next(s,:), 'fro' )  ), 0);
			R = R - G(:,s) * ( Y_next(s,:) - Y_prev(s,:) );
			Y_prev(s,:) = Y_next(s,:);
		end

		% if DEBUG == true
		% 	fprintf('delta = %g,', norm(Y_next - Y_prev, inf));
		% end

		% if norm(Y_next - Y_prev, inf) < tau
		% 	fprintf('breaked BCD, delta it = %d, MxNE_it = %d\n', i, k);
		% 	break;
		% end
		% Should compute a dual gap here and check for the convergence
% ------------------------------------------------------------------------------------ %
		eta  = dual_gap(M, G, Y_next, lambda, S, R);
% ------------------------------------------------------------------------------------ %
		

		if mod(i,5000) == 0 
				mu = mu / 2;
		end
		if eta < epsilon
			fprintf('breaked BCD, dual it = %d, MxNE_it = %d\n', i, k);
			break;
		elseif eta > 1e8
			fprintf('BCD ERROR: diverged!!\n');
			mu = mu /10;
			% break;
		end
	end
	fprintf('eta_bcd = %f, ', eta );
	Y = Y_next;
	iter = i;