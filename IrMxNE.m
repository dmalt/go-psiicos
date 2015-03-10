%% IrMxNE: linear regression with mixed norm regularization
%% function X = IrMxNE(M, G)
	
	% Initialization
	S = 100; 	% Number of sources
	T = 50; 	% Number of samples
	Sen = 11;	% Number of sensors

	X_prev = zeros(S, T);	% Matrix of signal sources; 
	X_next = zeros(S, T);	% Matrix of signal sources; 
	w = ones(S, 1); 	% Init weights vector
	W = eye(S); 		% Init weights matrix 

	lambda = 1.1; 	% Regularization parameter
	epsilon = 1e-6;	% Dual gap threshold
	tau = 1e-4;  	% Tolerance 
	I = 1000;		% Number of BCD iterations per one MxNE iteration
	K = 40;%30;		% Number of MxNE iterations

	tic;
	for k = 1:K
		X_prev = X_next;
		W = inv(diag(w));
		G = G_orig * W;
		l = zeros(S, 1);
		mu = ones(S, 1);%/1000;
		% for s = 1:S
		l = diag(G' * G);
		mu = 1. ./ ( 20. * l );
		% end

		Y_next = zeros(S, T);
		Y_prev = zeros(S, T);
		R = M;
		eta = epsilon * 2;

		for i = 1:I
			for s = 1:S
				Y_prev = Y_next;
				Y_next(s,:) = Y_prev(s,:) + mu(s) .* G(:,s)' * R;
				% if norm( Y_prev(s,:), 2 ) ~= 0
				Y_next(s,:) = Y_next(s,:)  - Y_prev(s,:) * mu(s) * lambda / (norm( Y_next(s,:), 2 )  );
				% else
					% fprintf('%d\n',s);
				% end
				R = R - G(:,s) * ( Y_next(s,:) - Y_prev(s,:) );

			end
			
			if norm(Y_next - Y_prev, inf) < tau
				fprintf('breaked BCD, it = %d, MxNE_it = %d\n', i, k);
				break;
			end
			% Should compute a dual gap here and check for the convergence 
		end
		X_next = W * Y_next;
		for s = 1:S
			w(s) = 1. / (  2 * sqrt(  norm( X_next(s,:), 2 )  )   );
		end
		if norm(X_next - X_prev, inf) < tau
			fprintf('breaked MxNE, it = %d\n', k);
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