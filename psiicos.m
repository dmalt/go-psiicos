%% IrMxNE: linear regression with mixed norm regularization
function X = IrMxNE(M, G)
	
	% Initialization
	S = 10 % Number of sources
	T = 10 % Number of samples

	X = zeros(S, T); % Matrix of signal sources; 
	w = ones(S, 1); % Init weights vector
	W = eye(S); % Init weights matrix 

	lambda = 1. % Regularization parameter
	epsilon = 1. % Dual gap threshold
	tau = 1e-6 % Tolerance 
	I = 100 % Number of BCD iterations per MxNE iteration
	K = 100 % Number of MxNE iterations

	for k = 1:K
		W = inv(diag(w));
		G = G * W;
		l = zeros(S,1);
		for s = 1:S
			fprintf("hi\n");
		end
	end

	X = 1;