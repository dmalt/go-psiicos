Init % Initialization procedure, where we set matrix dimensions

G = rand(Sen, S) ; % Gain matrix
G_orig = G;
Answer = zeros(S, T);
ix = randint(1,test_set_size, [1,S]);
Answer(ix,:) = rand(test_set_size, T) ;
M_signal = G * Answer; %rand(Sen, T); % Measurements. Should be passed from the outside
M = awgn(M_signal, SNR, 'measured');
% SNR = (norm(M_clear, 'fro') ^ 2 ) / (norm(M_noise,'fro') ^ 2)