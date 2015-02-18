% Initialization
%% IrMxNE: linear regression with mixed norm regularization
function X = IrMxNE(M, G)

	S = 10 % Number of sources
	T = 10 % Number of samples

	X = zeros(S, T); % Init our matrix of interest
	W = diag( ones(S * T, 1) ); % Init weights matrix
	w = ones(S * T, 1); % Weights vector

	lambda = 1. % Regularization parameter
	epsilon = 1. % Dual gap threshold
	tau = 1e-6 % Tolerance 
	I = 100 % Number of BCD iterations per MxNE iteration
	K = 100 % Number of MxNE iterations



	X = ;