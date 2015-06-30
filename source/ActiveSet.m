%% subset: returns an array of indices of sources which violate KKT conditions the most
function A = ActiveSet(G_small, R_, lambda, w, k, V)
	fprintf('Calculating active set...\n');
	fprintf('Start\n');
	startBCD = tic;
	R = [real(R_), imag(R_)];
	if k == 1		
		idx = @(x) x;
        S = length(w);
		save ../aux/G_small.txt G_small -ASCII
		save ../aux/R.txt R -ASCII
		unix('./rActSetViol');
		src_violations = load('../aux/V.txt');
		fprintf('Stop. ');
		src_violations = src_violations' - lambda;
	elseif k ~= 1
		idx = support(w);
		[dummy, S] = size(idx);
				src_violations(S) = 0; % Memory allocation. 
		for s = 1:S
			src_violations(s)  = norm(G_pair(idx(s),G_small, w(idx(s)))' * R, 'fro') - lambda;
		end	
	end
	[Dummy, A] =  sort(src_violations, 'descend');
	if V < S
		 A = A(1:V);
	end
	
	% src_violations(1:100)
	A = A(src_violations(A) > 0);
	% src_violations(A)
	A = idx(A);
	toc(startBCD);
	fprintf('Done\n');
	A = sort(A);
