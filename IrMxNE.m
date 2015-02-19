%% IrMxNE: linear regression with mixed norm regularization
%% function X = IrMxNE(M, G)
	
	% Initialization
	S = 50; 	% Number of sources
	T = 25; 	% Number of samples
	Sen = 10;	% Number of sensors

	X_prev = zeros(S, T);	% Matrix of signal sources; 
	X_next = zeros(S, T);	% Matrix of signal sources; 
	w = ones(S, 1); 	% Init weights vector
	W = eye(S); 		% Init weights matrix 

	lambda = 100.; 	% Regularization parameter
	epsilon = 1e-6;	% Dual gap threshold
	tau = 1e-6;  	% Tolerance 
	I = 100;		% Number of BCD iterations per MxNE iteration
	K = 500;		% Number of MxNE iterations

	G = rand(Sen, S)*100; % Gain matrix
	Answer = zeros(S, T);
	ix = randint(1,5, [1,S]);
	Answer(ix,:) = rand(5, T) * 100.;
	M = G * Answer; %rand(Sen, T); % Measurements. Should be passed from the outside

	for k = 1:K
		X_prev = X_next;
		W = inv(diag(w));
		G = G * W;
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
				Y_next(s,:) = Y_next(s,:) * (  1. - mu(s) * lambda / norm( Y_next(s,:), 2 )  );
				R = R - G(:,s) * ( Y_next(s,:) - Y_prev(s,:) );
			end
			% Should compute a dual gap here and check for the convergence 
		end
		X_next = W * Y_next;
		for s = 1:S
			w(s) = 1. / (   2. * sqrt(  norm( X_next(s,:) )  )   );
		end
		if norm(X_next - X_prev, inf) < tau
			fprintf('breaked');
			break;
		end
	end
	
			% Should compute dual gap here and check for convergence 
	norm(X_next, 'fro')
	norm(Answer, 'fro')
	norm(X_next - Answer, 'fro')