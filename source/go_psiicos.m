
function [X,Aidx] = go_psiicos(CT, CT2, G_small, lambda, ActSetChunk, Subj)
% go_psiicos: compute connectivity on sources given cross-spectrum on sensors 
%             and forward model matrix solving a global optimization problem. 
% FORMAT:
%   [X,Aidx] = go_psiicos(CT, CT2, G_small, lambda, ActSetChunk, Subj)
% INPUTS:
%   CT                 - {Nsensors^2 x Time} cross-spectrum on sensors
%   CT2                - {Nsensors^2 x Time} cross-spectrum we want project from
%   G_small            - {Nsensors x Nsources} forward model matrix in physical space
%   lambda             - regularization parameter (default = 0.05)
%   ActSetChunk        - chunk of sources to be added to the active set on each iteration
%   Subj               - name of a subject. If simulations, assign 'simul'.
% OUTPUTS:
%   X                  - {N_act_sources x Time} recovered timecources on sources
%   Aidx               - indices of recovered connections 
% __________________________________________________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
 
    % Initialization
    DEBUG = 0;
% -------------------------------------------------------------------------- %
    % M_abs = abs(M);
    setenv('OPENBLAS_NUM_THREADS', '1');
    setenv('OMP_NUM_THREADS', '4');
    if nargin < 5 
        CT2 = [];
    end
    %--- Reduce time timecourse lengthes using pca ----------------------------%
    M_full_complex  = ProjOut(CT, CT2, G_small) ;
    M_norm = M_full_complex / norm(M_full_complex);
    ncomp = 8;
    Re_M_norm = real(M_norm);
    Im_M_norm = imag(M_norm);
    M_cat = [Re_M_norm, Im_M_norm];
    [Mu Ms Mv] = svd(M_cat);
    M_cat_reduced = M_cat * Mv(:,1:ncomp); % M_cat_reduced is real matrix with timecourses of double length
    Tx2_full = size(M_cat,2);
    T_full = Tx2_full / 2;
    T = size(M_cat_reduced,2); % Nch - number of channels, T - number of time samples
    % T = floor(Tx2 / 2);
    M = M_cat_reduced;%(:,1:T) + j * M_cat_reduced(:,T+1 : Tx2);
    clear M_full_complex M_norm Re_M_norm Im_M_norm M_cat M_cat_reduced;
    % -------------------------------------------------------------------------%
    [Nch, Nsrc] = size(G_small);
    Nsrc_pairs = Nsrc ^ 2; 
    Nsites = Nsrc / 2; % One site contains two dipoles
    Nsite_pairs = Nsites ^ 2; % Pairs of sites

    % Normalize generating matix %
    for i = 1:Nsites
         range_i = i*2-1:i*2;
         G_small(:,range_i(1)) = G_small(:,range_i(1))/norm(G_small(:,range_i(1)));
         G_small(:,range_i(2)) = G_small(:,range_i(2))/norm(G_small(:,range_i(2)));
    end;

    A_reduced = [];
    X_prev_active = []; % Matrix of signal sources; 
    X_next_active = []; % Matrix of signal sources; 
    suppX_next = [];
    suppX_prev = [];
    G_a = zeros(Nch* Nch , 4);
    X_a = zeros(4, T);
    w = ones(Nsite_pairs, 1);   % Init weights vector

    % lambda = 0.08;    % Regularization parameter
    epsilon = 1e-5;     % Dual gap threshold
    eta = 2;            % Primal-dual gap  
    tau = 1e-4;         % Tolerance 
    K = 50;%30;         % Number of MxNE iterations

% --------------------------------------------------------------------------- %
    timerOn = tic;
    for k = 1:K
% ------------------------------------------------------------------------------- %
        % ActSetChunk = 500;
        Res = M;
        A = ActiveSet(G_small, Res, lambda, w, k, ActSetChunk);
        if isempty(A)
            break;
        end
        fprintf('Calculating max l(s)...\n');
        l(Nsite_pairs) = 0;
    % To accelerate calculation on the second and subsequent iterations making use of the reduced structurre of G when k~=1
        if k == 1 
            S = Nsite_pairs; 
            % if matlabpool('size') == 0
                matlabpool open;
            % end
            parfor s = 1:S
                G_s = G_pair(  s, G_small, w(s)  );
                l(s) = norm(G_s' * G_s, 'fro');
            end
            % if matlabpool('size') > 0
                matlabpool close;
            % end
        elseif k ~= 1
            idx = support(w);
            [dummy, S] = size(idx);  
            for s = 1:S
                G_s = G_pair(  idx(s), G_small, w( idx(s) )  );
                l(s) = norm(G_s' * G_s, 'fro');
            end
        end
        mu = 1 / max(l);
         % mu = 0.499257
        fprintf('Done.\n'); 
        fprintf('mu = %f\n', mu );
        % mu(support(l)) = 1. ./  (5*l(support(l)));    
 %  ------------------------------------------------------------------------------ %
        % mu = ones(Nsrc_pairs, 1) / 1000;
        % X = zeros(Nsrc_pairs, T);
        
        % if k == 1
        %   while isempty(A) 
        %       lambda = lambda / 2;
        %       fprintf('lambda reduced to %f\n', lambda);
        %       A = ActiveSet(G_small, Res, lambda, w, k, ActSetChunk);
        %   end
        % elseif k ~= 1
        %   if isempty(A)
        %       fprintf('A is empty');
        %       break;
        %   end
        % end
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
            [X_a, bcd_iter] = BCD__(G_a, X_a, M, lambda, epsilon, mu);
            % X = zeros(Nsrc_pairs,T);
            A_reduced = A(1, supp_d(X_a)); 
            if k ~= 1
                disp([(mod(A,Nsites))', ((A - mod(A,Nsites)) / Nsites+ 1)']);
            end
            X_a_reduced = X_a(ind4(supp_d(X_a)),:);
            Res = M - G_a * X_a;
            % X(A_reduced,:) = X_a_reduced;
            fprintf('BCD iter = %d, eta = %f, Active set size = %d\n', bcd_iter, eta, sizeA);
        
            A_penalized = ActiveSet(G_small, Res, lambda, w, k, ActSetChunk);
            eta = dual_gap_big([real(M), imag(M)], G_small, [real(X_a), imag(X_a)], lambda, sizeA, [real(Res), imag(Res)]);
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
            save ../output/A_reduced.mat A_reduced
        end
% ------------------------------------------------------------------------------------------ %
        suppX_next = A_reduced;
        w_ = bsxfun(@times,w(suppX_next,1)',  ones(4,1)); % Making 4 copies of weight to calc new X
        w_ = w_(:); 
        X_next_active = diag(w_) * X_a_reduced; 
        % for s = 1:Nsrc_pairs
        %   w(s) =   2 * sqrt(  norm( X_next(s,:), 2 )  ) ;
        % end
        w = zeros(Nsite_pairs,1);
        Sw = length(suppX_next);
        for sw = 1:Sw
            w(suppX_next(sw),1) = 2 * sqrt(  norm( X_next_active(ind4(sw),:), 'fro' )  ); 
        end

        % --- Calculate norm of difference between solutions on current and previous step ---%
        suppUnion = sort(union(suppX_next, suppX_prev));
        sizeSuppUnion = length(suppUnion); 
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

    N1 =  norm(M - G_a * X_a, 'fro');
    actSites = size(X_next_active,1) / 4;
    N2 = 0;
    range = 1:4;
    for s = 1:actSites
        N2 = N2 + norm(X_a(range,:), 'fro');
        range = range + 4; 
    end
    if matlabpool('size') > 0
        matlabpool close;
    end
    exec_time = toc(timerOn);
    lambda_str = num2str(lambda);
    % ReX_next_active = real(X_next_active);
    % ImX_next_active = imag(X_next_active);
    % X_next_active_cat = [ReX_next_active, ImX_next_active];
    X_rec_cat = X_next_active * Mv(:,1:ncomp)';
    X = X_rec_cat(:,1:T_full) + j * X_rec_cat(:,T_full+1:Tx2_full); 
    Aidx = A_reduced;
    save (  strcat('../output/Output_', Subj,'_', lambda_str,'_', date, '.mat'), 'Aidx','X', 'N1', 'N2', 'lambda', 'CT', 'CT2','G_small', 'epsilon', 'tau', 'ncomp', 'ActSetChunk', 'exec_time', 'i', 'k');
end
