%% IrMxNE: linear regression with mixed norm regularization
%% function X = IrMxNE(M, G)
	
	% Initialization
	S = 100; 	% Number of sources
	T = 50; 	% Number of samples
	Sen = 25;	% Number of sensors

	X_prev = zeros(S, T);	% Matrix of signal sources; 
	X_next = zeros(S, T);	% Matrix of signal sources; 
	w = ones(S, 1); 	% Init weights vector
	W = eye(S); 		% Init weights matrix 

	lambda = 2; 	% Regularization parameter
	epsilon = 1e-6;	% Dual gap threshold
	tau = 1e-6;  	% Tolerance 
	I = 100;		% Number of BCD iterations per MxNE iteration
	K = 100;%30;		% Number of MxNE iterations

	for k = 1:K
		X_prev = X_next;
		W = inv(diag(w));
		G = G_orig * W;
		l = zeros(S, 1);
		mu = zeros(S, 1);
		for s = 1:S
			l(s) = G(:, s)' * G(:, s);
			mu(s) = 1. / ( 2. * l(s) );
		end

		Y_next = zeros(S, T);
		Y_prev = zeros(S, T);
		R = M;
		eta = epsilon * 2;

		for i = 1:I
			for s = 1:S
				Y_prev = Y_next;
				Y_next(s,:) = Y_prev(s,:) + mu(s) .* G(:,s)' * R;
				if norm( Y_prev(s,:), 2 ) ~= 0
					Y_next(s,:) = Y_next(s,:)  - Y_prev(s,:) * mu(s) * lambda / (norm( Y_prev(s,:), 2 )  );
				% else
					% fprintf('%d\n',s);
				end
				R = R - G(:,s) * ( Y_next(s,:) - Y_prev(s,:) );
			end
			% Should compute a dual gap here and check for the convergence 
		end
		X_next = W * Y_next;
		for s = 1:S
			w(s) = 1. / (  2 * sqrt(  norm( X_next(s,:), 2 )  )   );
		end
		if norm(X_next - X_prev, inf) < tau
			fprintf('breaked\n');
			break;
		end
	end
	
			% Should compute dual gap here and check for convergence 
	norm(X_next, 'fro')
	norm(Answer, 'fro')
	norm(X_next - Answer, 'fro')
	norm(M - G_orig * X_next, 'fro')
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