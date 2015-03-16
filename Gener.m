Init

G = rand(Sen, S) ; % Gain matrix
G_orig = G;
Answer = zeros(S, T);
ix = randint(1,test_set_size, [1,S]);
Answer(ix,:) = rand(test_set_size, T) ;
M = G * Answer; %rand(Sen, T); % Measurements. Should be passed from the outside