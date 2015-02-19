%% IrMxNE: linear regression with mixed norm regularization
%function X = IrMxNE(M, G)
	
	% Initialization
	S = 10; 	% Number of sources
	T = 10; 	% Number of samples
	Sen = 2;	% Number of sensors

	X = zeros(S, T);	% Matrix of signal sources; 
	w = ones(S, 1); 	% Init weights vector
	W = eye(S); 		% Init weights matrix 

	lambda = 1.; 	% Regularization parameter
	epsilon = 1e-6;	% Dual gap threshold
	tau = 1e-6;  	% Tolerance 
	I = 100;		% Number of BCD iterations per MxNE iteration
	K = 100;		% Number of MxNE iterations

	M = zeros(Sen, T) % Measurement. M should be passed from outside
	G = zeros(Sen, S) % Gain matrix

	for k = 1:K
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
				Y_next(s,:) = Y_prev(s,:) + mu(s) .* G(:,s)' * R;
				Y_next(s,:) = Y_next(s,:) * (  1 - mu(s) * lambda / norm( Y_next(s,:), 2 )  );
				R = R - G(:,s) * ( Y_next(s,:) - Y_prev(s,:) );
			end
		end

	end
	X = 1;