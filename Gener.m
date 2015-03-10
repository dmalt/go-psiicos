S = 200; 	% Number of sources
T = 500; 	% Number of samples
Sen = 20;	% Number of sensors

G = rand(Sen, S) ; % Gain matrix
G_orig = G;
Answer = zeros(S, T);
ix = randint(1,5, [1,S]);
Answer(ix,:) = rand(5, T) ;
M = G * Answer; %rand(Sen, T); % Measurements. Should be passed from the outside