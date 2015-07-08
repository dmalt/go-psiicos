%% subset: returns an array of indices of sources which violate KKT conditions the most
function A = ActiveSet(G_small, R_, lambda, w, k, V)
	fprintf('Calculating active set...\n');
	fprintf('Start\n');
	tic;
	R = [real(R_), imag(R_)];
	save ../aux/G_small.txt G_small -ASCII
	save ../aux/R.txt R -ASCII
	save ../aux/w.txt w -ASCII
	unix('./ActSetViol');
	src_violations = load('../aux/V.txt');
	fprintf('Stop. ');
	toc;
	src_violations = src_violations' - lambda;
	[Dummy, A] =  sort(src_violations, 'descend');
	 % if 2*V < S
		 A = A(1:V);
	 % end
	
	% src_violations(1:100)
	A = A(src_violations(A) > 0);
	% src_violations(A)
	% A = idx(A);
	fprintf('Done\n');
	A = sort(A);
