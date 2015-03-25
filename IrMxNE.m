%% IrMxNE: linear regression with mixed norm regularization; M - measurements, G_small - forward model matrix in physical space
function X = IrMxNE(M, G_small)
	DEBUG = true;

	% Initialization
% -------------------------------------------------------------------------- %
	[Nch , T] = size(M); % Nch - number of channels, T - number of time samples
	[Nch_small, Nsrc_small] = size(G_small);
	Nsrc = Nsrc_small ^ 2;

	X_prev = zeros(Nsrc, T);	% Matrix of signal sources; 
	X_next = zeros(Nsrc, T);	% Matrix of signal sources; 
	w = ones(Nsrc, 1); 	% Init weights vector

	lambda = 80.; 		% Regularization parameter
	epsilon = 1e-5;		% Dual gap threshold
	eta = epsilon * 2;	% Primal-dual gap  
	tau = 1e-4;  		% Tolerance 
	K = 50;%30;			% Number of MxNE iterations
% --------------------------------------------------------------------------- %
	tic;
	for k = 1:K
		X_prev = X_next;
		l = zeros(Nsrc, 1);
		% figure; image(G*100);
		% l = diag(G' * G);				% !!!! FIX THIS !!!! 
		mu = ones(size(l))/max(l);		% !!!! FIX THIS !!!!
		% mu(support(l)) = 1. ./  (5*l(support(l)));
 %  ------------------------------------------------------------------------------ %
			% mu = ones(Nsrc, 1) / 1000;
			% X = zeros(Nsrc, T);
			R = M;
			A = ActiveSet(G_small, R, lambda);
			A_reduced = [];
			X_a = zeros(A, T);
			X_a_reduced = X_a;
			for i = 1:100
				G_a = CalcG(A, G_small, w);
				[dummy, sizeA] = size(A);
				X_a = zeros(A, T);
				[dummy, nonzero_idx] = intersect(A, A_reduced);
				X_a(nonzero_idx, :) =  X_a_reduced;
				[dummy, size_A] = size(A);
				[X_a, bcd_iter] = BCD(size_A, T, G_a, X_a, M, lambda, epsilon, k, mu, tau, DEBUG );
				% X = zeros(Nsrc,T);
				A_reduced = A(1, support(X_a));
				X_a_reduced = X_a(support(X_a),:)
				R = M - G_a * X_a;
				% eta = dual_gap(M, G, X, lambda, Nsrc, R);
				fprintf('bcd iterations = %d, eta = %f, size_A = %d\n', bcd_iter, eta, size_A);
			
				A_penalized = ActiveSet(G_small, R, lambda);
				A_next = sort(union(A_reduced, A_penalized));
				isthesame = isempty(setxor(A, A_next));
				if isthesame
					fprintf('isempty = %d, eta = %f, ', isthesame, eta );
					break;
				end
				A = A_next;
			end
% -------------------------------------------------------------------------------- %
			X_next = W * X; % !!!!!!!!!!!!!!!!!!!!!!! FIX THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			% for s = 1:Nsrc
			% 	w(s) =   2 * sqrt(  norm( X_next(s,:), 2 )  ) ;
			% end
			w = 2 * sqrt(diag(X_next * X_next'));

			fprintf('delta_1 = %g\n', norm(X_next - X_prev, inf));
			
		if norm(X_next - X_prev, inf) < tau
			fprintf('breaked from irMxNE, it = %d\n', k);
			break;
		end
	end
	elapsed = toc;
	fprintf('TIC TOC: %g\n', elapsed);
			% Should compute dual gap here and check for convergence 
	Result_norm = norm(X_next, 'fro');
	Answer_norm = norm(Answer, 'fro');
	Error_norm = norm(X_next - Answer, 'fro');
	Residual_norm = norm(M - G_orig * X_next, 'fro');
	fprintf('Result_norm: %g\n', Result_norm);
	fprintf('Answer_norm: %g\n', Answer_norm);
	fprintf('Error_norm: %g\n', Error_norm);
	fprintf('Residual_norm: %g\n', Residual_norm);
	F1 = 2 * length( intersect(support(X_next), support(Answer) ) ) / (length(support(X_next)) + length(support(Answer)));
	fprintf('F1 =  %f\n', F1 );
	if isempty(setxor(support(X_next),ix))
		fprintf('Hooray, we`ve found all the sources\n');
	end
	% [I,J] = size(X_next);
	% for i = 1:I
	% 	for j = 1:J
	% 		if X_next(i,j) < 0.11
	% 			X_next(i,j) = 0;
	% 		end
	% 	end
	% end
	% norm(X_next - Answer, 'fro')
	% norm(M - G_orig * X_next, 'fro')