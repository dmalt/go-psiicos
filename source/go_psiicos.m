%% IrMxNE: linear regression with mixed norm regularization; M - measurements, G_small - forward model matrix in physical space
function [X,Aidx] = IrMxNE(G_small, CT, CT2)
	% Initialization
    DEBUG = 0;
% -------------------------------------------------------------------------- %
	% M_norm = abs(M);
	if nargin < 3 
		CT2 = [];
	end
	M  = ProjOut(CT, CT2, G_small) ;
	M_norm = M / norm(M);
	ncomp = 5;
	[Mu Ms Mv] = svd(M_norm);
	M_svd_full = M_norm * Mv(:,1:ncomp);
	
	[dummy , T] = size(M_svd_full); % Nch - number of channels, T - number of time samples
	[Nch, Nsrc] = size(G_small);
	Nsrc_pairs = Nsrc ^ 2; 
	Nsites = Nsrc / 2; % One site contains two dipoles
	Nsite_pairs = Nsites ^ 2; % Pairs of sites

	IND((Nch ^ 2 + Nch) / 2) = 0; 
	s = 1;
	for k = 1:Nch
    	for l = k:Nch
        	IND(s) = Nch * (k - 1) + l;
        	s = s + 1;
    	end
	end
	M_svd = M_svd_full(IND,:);
	% Normalize generating matix %
	for i = 1:Nsites
   		 range_i = i*2-1:i*2;
   		 G_small(:,range_i(1)) = G_small(:,range_i(1))/norm(G_small(:,range_i(1)));
   		 G_small(:,range_i(2)) = G_small(:,range_i(2))/norm(G_small(:,range_i(2)));
	end;

	X_prev_active = [];	% Matrix of signal sources; 
	X_next_active = [];	% Matrix of signal sources; 
	suppX_next = [];
	suppX_prev = [];
	w = ones(Nsite_pairs, 1); 	% Init weights vector

	lambda = 0.005;		% Regularization parameter

	epsilon = 1e-5;		% Dual gap threshold
	eta = 2;			% Primal-dual gap  
	tau = 1e-4;  		% Tolerance 
	K = 50;%30;			% Number of MxNE iterations

% --------------------------------------------------------------------------- %
	timerOn = tic;
	for k = 1:K
% ------------------------------------------------------------------------------- %
		ActSetChunk = 25;
		fprintf('Calculating max l(s)...\n');
		l(Nsite_pairs) = 0;
	% To accelerate calculation on the second and subsequent iterations making use of the reduced structurre of G when k~=1
		if k == 1 
			S = Nsite_pairs; 
			idx = @(x) x;
			matlabpool('open', 4);
			parfor s = 1:S
				G_s = G_pair(  idx(s), G_small, w( idx(s) )  );
				l(s) = norm(G_s' * G_s, 'fro');
			end
			matlabpool close;
		elseif k ~= 1
			idx = support(w);
			[dummy, S] = size(idx);  
			for s = 1:S
				G_s = G_pair(  idx(s), G_small, w( idx(s) )  );
				l(s) = norm(G_s' * G_s, 'fro');
			end
		end
		
		mu = 1 / max(l);

		fprintf('Done.\n');	
		fprintf('mu = %f\n', mu );
		% mu(support(l)) = 1. ./  (5*l(support(l))); 	
 %  ------------------------------------------------------------------------------ %
		% mu = ones(Nsrc_pairs, 1) / 1000;
		% X = zeros(Nsrc_pairs, T);
		Res = M_svd;
		A = ActiveSet(G_small, Res, lambda, w, k, ActSetChunk);
		if k == 1
			while isempty(A) 
				lambda = lambda / 2;
				fprintf('lambda reduced to %f\n', lambda);
				A = ActiveSet(G_small, Res, lambda, w, k, ActSetChunk);
			end
		elseif k ~= 1
			if isempty(A)
				fprintf('A is empty');
				break;
			end
		end
		A_reduced = [];
		X_a_reduced = [];
		for i = 1:100
			[dummy, sizeA] = size(A);
			G_a = CalcG(A, G_small, w);
			X_a = zeros(sizeA * 4, T);
			[dummy, nonzero_idx] = intersect(A, A_reduced);
			if ~isempty(nonzero_idx)
				X_a(ind4(nonzero_idx), :) =  X_a_reduced;	
			end
			[X_a, bcd_iter] = BCD(G_a, X_a, M_svd, lambda, epsilon, mu);
			% X = zeros(Nsrc_pairs,T);
			A_reduced = A(1, supp_d(X_a)); 
			Pairs = [(mod(A,Nsites))', ((A - mod(A,Nsites)) / Nsites+ 1)']
			X_a_reduced = X_a(ind4(supp_d(X_a)),:);
			Res = M_svd - G_a * X_a;
			% X(A_reduced,:) = X_a_reduced;
			fprintf('BCD iter = %d, eta = %f, Active set size = %d\n', bcd_iter, eta, sizeA);
		
			A_penalized = ActiveSet(G_small, Res, lambda, w, k, ActSetChunk);
			eta = dual_gap_big([real(M_svd), imag(M_svd)], G_small, [real(X_a), imag(X_a)], lambda, sizeA, [real(Res), imag(Res)]);
			fprintf('eta = %f, ', eta );

			A_next = sort(union(A_reduced, A_penalized));
			isthesame = isempty(setxor(A, A_next));
			if isthesame
				 % ActSetChunk = ActSetChunk * 2;
			% end
			% if eta < epsilon
				fprintf('isempty = %d, eta = %f, ', isthesame, eta );
				break;
			end
			A = A_next;
			save ../output/A_reduced.mat A_reduced;
		end
% ------------------------------------------------------------------------------------------ %
		suppX_next = A_reduced;
		w_ = bsxfun(@times,w(suppX_next,1)',  ones(4,1)); % Making 4 copies of weight to calc new X
		w_ = w_(:); 
		X_next_active = diag(w_) * X_a_reduced; 
		% for s = 1:Nsrc_pairs
		% 	w(s) =   2 * sqrt(  norm( X_next(s,:), 2 )  ) ;
		% end
		w = zeros(Nsite_pairs,1);
		Sw = length(suppX_next);
		for sw = 1:Sw
			w(suppX_next(sw),1) = 2 * sqrt(  norm( X_next_active(ind4(sw),:), 'fro' )  ); 
		end

		% --- Calculate norm of difference between solutions on current and previous step ---%
		suppUnion = sort(union(suppX_next, suppX_prev));
		[dummy, sizeSuppUnion] = size(suppUnion); 
		X_next_exp = zeros(sizeSuppUnion * 4,T);
		X_prev_exp = zeros(sizeSuppUnion * 4,T);
		[dummy, nonzero_idx_next] = intersect(suppUnion, suppX_next);
		[dummy, nonzero_idx_prev] = intersect(suppUnion, suppX_prev);
		if ~isempty(nonzero_idx_next)
			X_next_exp(ind4(nonzero_idx_next),:) = X_next_active;
		end
		if ~isempty(nonzero_idx_prev)
			X_prev_exp(ind4(nonzero_idx_prev),:) = X_prev_active;
		end
		fprintf('delta_1 = %g\n', norm(X_next_exp - X_prev_exp, inf));
	
		if norm(X_next_exp - X_prev_exp, inf) < tau
			fprintf('breaked from irMxNE, it = %d\n', k);
			break;
		end
		% ---------------------------------------------------------------------------------- %
		X_prev_active = X_next_active;
		suppX_prev = suppX_next;
	end

	toc(timerOn);
	lambda_str = num2str(lambda);
	save ( strcat( strcat('../output/Output_', lambda_str), '.mat'), 'A_reduced','X_next_active');
    X = X_next_active;
    Aidx = A_reduced;
