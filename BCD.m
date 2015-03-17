function [Y, iter]  = BCD(S, T, G, Y_prev, M, lambda, epsilon, k, mu, tau, DEBUG)
	I = 10000;			% Number of BCD iterations per one MxNE iteration
	Y_next = zeros(S, T);
	R = M - G * Y_prev;

	for i = 1:I
		for s = 1:S
			
			Y_next(s,:) = Y_prev(s,:) + mu(s) .* G(:,s)' *	(M - G * Y_prev);
			% if norm( Y_prev(s,:), 2 ) ~= 0
			Y_next(s,:) = Y_next(s,:) * max(1 -  mu(s) * lambda / (norm( Y_next(s,:), 'fro' )  ), 0);
			% else
				% fprintf('%d\n',s);
			% end
			
			R = R - G * ( Y_next - Y_prev );%G(:,s) * ( Y_next(s,:) - Y_prev(s,:) );
			Y_prev = Y_next;
		end

		if DEBUG == true
			fprintf('delta = %g,', norm(Y_next - Y_prev, inf));
		end

		% if norm(Y_next - Y_prev, inf) < tau
		% 	fprintf('breaked BCD, delta it = %d, MxNE_it = %d\n', i, k);
		% 	break;
		% end
		% Should compute a dual gap here and check for the convergence
% ------------------------------------------------------------------------------------ %
		[eta, Fp, Fd]  = dual_gap(M, G, Y_next, lambda, S, R);
% ------------------------------------------------------------------------------------ %
		if DEBUG == true
			fprintf('bcd iter = %d, k = %d, eta = %f, Fp = %f, Fd = %f\n', i, k, eta, Fp, Fd);
		end


		if eta < epsilon
			fprintf('breaked BCD, dual it = %d, MxNE_it = %d\n', i, k);
			break;
		end
	end
	fprintf('eta_bcd = %f, ', eta );
	Y = Y_next;
	iter = i;