%% subset: returns an array of indices of sources which violate KKT conditions the most
function A = ActiveSet(G_small, R, lambda, w, k)
	fprintf('Calculating active set...\n');
	V = 200; % Number of sources in a subset
	% [ans, S] = size(G_small);
	% if k == 1
	% 	S = 2 * S ^ 2; 
	% 	idx = @(x) x;
	% elseif k ~= 1
	% 	idx = support(w);
	% 	[dummy, S] = size(idx);
	% end


	% matlabpool('open', 4);
	% src_violations(S) = 0; % Memory allocation. 
	% parfor s = 1:S
	% 	src_violations(s)  = norm(columnG_fast(idx(s),G_small, w)' * R) - lambda;
	% end	
	% matlabpool close;
	fprintf('Start\n');
	tic;
	save G_small.txt G_small -ASCII
	save R.txt R -ASCII
	save w.txt w -ASCII
	unix('./load_matr');
	src_violations = load('V.txt');
	fprintf('Stop. ');
	toc;
	src_violations = src_violations' - lambda;
	[Dummy, A] =  sort(src_violations, 'descend');
	% if V < S
		A = A(1:V);
	% end
	A = A(src_violations(A) > 0);
	% src_violations(A)
	% A = idx(A);
	fprintf('Done\n');
	A = sort(A);
