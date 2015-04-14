%% IrMxNE: linear regression with mixed norm regularization; M - measurements, G_small - forward model matrix in physical space
% function X = IrMxNE(M, G_small)
	DEBUG = true;

	% Initialization

% -------------------------------------------------------------------------- %
	M_real = [real(M); imag(M)];
	G_small = G2dU;
	[Nch , T] = size(M_real); % Nch - number of channels, T - number of time samples
	[Nch_small, Nsrc_small] = size(G_small);
	Nsrc = 2 * Nsrc_small ^ 2; % Because we want to get real instead of dealing with a complex space


	X_prev_active = [];	% Matrix of signal sources; 
	X_next_active = [];	% Matrix of signal sources; 
	suppX_next = [];
	suppX_prev = [];
	w = ones(Nsrc, 1); 	% Init weights vector

	lambda = 250.; 		% Regularization parameter
	epsilon = 1e-5;		% Dual gap threshold
	eta = 2;	% Primal-dual gap  
	tau = 1e-4;  		% Tolerance 
	K = 50;%30;			% Number of MxNE iterations
% --------------------------------------------------------------------------- %
	tic;
	for k = 1:K
		% X_prev = X_next;
		% l = zeros(Nsrc, 1);
		% figure; image(G*100);
		fprintf('Calculating max l(s)...\n');
		l(Nsrc) = 0;
		if k == 1
			S = Nsrc; 
			idx = @(x) x;
		elseif k ~= 1
			idx = support(w);
			[dummy, S] = size(idx);
		end
		matlabpool('open', 4);
		parfor s = 1:S
			l(s) = sum(columnG_fast(idx(s), G_small, w).^2);
		end
		matlabpool close;
		mu = 1 / max(l);
		fprintf('Done.\n');	
		fprintf('mu = %f\n', mu );
		% mu(support(l)) = 1. ./  (5*l(support(l))); 	
 %  ------------------------------------------------------------------------------ %
		% mu = ones(Nsrc, 1) / 1000;
		% X = zeros(Nsrc, T);
		R = M_real;
		
		A = ActiveSet(G_small, R, lambda, w, k);
		A_reduced = [];
		X_a_reduced = [];
		for i = 1:100
			G_a = CalcG(A, G_small, w);
			[dummy, sizeA] = size(A); 
			X_a = zeros(sizeA, T);
			[dummy, nonzero_idx] = intersect(A, A_reduced);
			if ~isempty(nonzero_idx)
				X_a(nonzero_idx, :) =  X_a_reduced;	
			end
			[X_a, bcd_iter] = BCD(sizeA, T, G_a, X_a, M_real, lambda, epsilon, k, mu, tau, DEBUG );
			% X = zeros(Nsrc,T);
			A_reduced = A(1, support(X_a));
			X_a_reduced = X_a(support(X_a),:);
			R = M_real - G_a * X_a;
			% eta = dual_gap(M_real, G, X, lambda, Nsrc, R);
			fprintf('BCD iter = %d, eta = %f, Active set size = %d\n', bcd_iter, eta, sizeA);
		
			A_penalized = ActiveSet(G_small, R, lambda, w, k);
			A_next = sort(union(A_reduced, A_penalized));
			isthesame = isempty(setxor(A, A_next));
			if isthesame
				fprintf('isempty = %d, eta = %f, ', isthesame, eta );
				break;
			end
			A = A_next;
		end
% ------------------------------------------------------------------------------------------ %
		suppX_next = A_reduced;
		X_next_active = diag(w(suppX_next,1)) * X_a_reduced; 
		% for s = 1:Nsrc
		% 	w(s) =   2 * sqrt(  norm( X_next(s,:), 2 )  ) ;
		% end
		w = zeros(Nsrc,1);
		w(suppX_next,1) = 2 * sqrt(diag(X_next_active * X_next_active')); 
		% --- Calculate norm of difference between solutions on current and previous step ---%
		suppUnion = sort(union(suppX_next, suppX_prev));
		[dummy, sizeSuppUnion] = size(suppUnion); 
		X_next_exp = zeros(sizeSuppUnion,T);
		X_prev_exp = zeros(sizeSuppUnion,T);
		[dummy, nonzero_idx_next] = intersect(suppUnion, suppX_next);
		[dummy, nonzero_idx_prev] = intersect(suppUnion, suppX_prev);
		if ~isempty(nonzero_idx_next)
			X_next_exp(nonzero_idx_next,:) = X_next_active;
		end
		if ~isempty(nonzero_idx_prev)
			X_prev_exp(nonzero_idx_prev,:) = X_prev_active;
		end
		fprintf('delta_1 = %g\n', norm(X_next_exp - X_prev_exp, inf));
	
		if norm(X_next_exp - X_prev_exp, inf) < tau
			fprintf('breaked from irMxNE, it = %d\n', k);
			break;
		end
		% ---------------------------------------------------------------------------------- %
		X_prev_active = X_next_active;
		suppX_prev = suppX_next;
	end
	elapsed = toc;
	fprintf('TIC TOC: %g\n', elapsed);
			% Should compute dual gap here and check for convergence 
	% [I,J] = size(X_next);
	% for i = 1:I
	% 	for j = 1:J
	% 		if X_next(i,j) < 0.11
	% 			X_next(i,j) = 0;
	% 		end
	% 	end
	% end
	% norm(X_next - Answer, 'fro')
	% norm(M_real - G_orig * X_next, 'fro')
