%% IrMxNE: linear regression with mixed norm regularization
%% function X = IrMxNE(M, G)
	DEBUG = true;

	% Initialization
% -------------------------------------------------------------------------- %
	Init

	X_prev = zeros(S, T);	% Matrix of signal sources; 
	X_next = zeros(S, T);	% Matrix of signal sources; 
	w = ones(S, 1); 	% Init weights vector
	W = eye(S); 		% Init weights matrix 

	lambda = 800.; 		% Regularization parameter
	epsilon = 1e-5;		% Dual gap threshold
	eta = epsilon * 2;	% Primal-dual gap  
	tau = 1e-4;  		% Tolerance 
	K = 50;%30;			% Number of MxNE iterations
% --------------------------------------------------------------------------- %
	tic;
	for k = 1:K
		X_prev = X_next;
		W = diag(w);
		G = G_orig * W;
		l(S) = 0;
		% figure; image(G*100);
		fprintf('Calculating l...\n');
		l = diag(G' * G);
		mu = 1/max(l);
		fprintf('Done.\n');
		% mu(support(l)) = 1. ./  (5*l(support(l)));
 %  ------------------------------------------------------------------------------ %
			% mu = ones(S, 1) / 1000;
			X = zeros(S, T);
			R = M;
			A = subset(G, R, lambda, S);
			
			for i = 1:100
				G_a = G(:,A);
				X_a = X(A,:);
				[dummy, size_A] = size(A);
				[X_a, bcd_iter] = BCD(size_A, T, G_a, X_a, M, lambda, epsilon, k, mu, tau, DEBUG );
				X = zeros(S,T);
				X(A,:) = X_a;
				R = M - G * X;
				eta = dual_gap(M, G, X, lambda, S, R);
				fprintf('bcd iterations = %d, eta = %f, size_A = %d\n', bcd_iter, eta, size_A);
			
				A_ = subset(G, R, lambda, S);
				A_next = union(support(X), A_);
				[dummy, size_A] = size(A);
				[dummy, size_A_next] = size(A_next);
				isthesame = isempty(setxor(A,A_next));
				if eta < epsilon || isthesame
					fprintf('isempty = %d, eta = %f, ', isthesame, eta );
					break;
				end
				A = A_next;
			end
% -------------------------------------------------------------------------------- %
			X_next = W * X;
			for s = 1:S
				w(s) =   2 * sqrt(  norm( X_next(s,:), 2 )  ) ;
			end

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