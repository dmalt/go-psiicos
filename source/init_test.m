%% init_test: Generates simple test data for algorythm validation
function [M, G2DU, Answer, ans_idx] = init_test(Src, Sen, T, SNR, test_set_size)

% Initialization procedure, where we set matrix dimensions
% Src = 50; 	% Number of sources
% T = 10; 	% Number of samples
% Sen = 20;	% Number of sensors
% SNR = 2;	% Signal to noise ratio
% test_set_size = 10;
G2DU = rand(Sen, Src) ; % Gain matrix
G_orig = G2DU;
Answer = zeros(Src^2, T);
G_big = [];
for i = 1:Src
	for j = 1:Src
		matGbig = G2DU(:,i) * G2DU(:,j)';
		colGbig = matGbig(:); 
		G_big = [G_big, colGbig];
	end
end
ans_idx = randint(1, test_set_size, [1, Src ^ 2]);
Answer(ans_idx,:) = ones(test_set_size, T) * 2;
M_signal = G_big * Answer; %rand(Sen, T); % Measurements. Should be passed from the outside
M = awgn(M_signal, SNR, 'measured');
% SNR = (norm(M_clear, 'fro') ^ 2 ) / (norm(M_noise,'fro') ^ 2)
