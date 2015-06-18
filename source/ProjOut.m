function [ Cp ] = ProjOut(CT, CT2, G)

% if(nargin < 3)
    Rnk  = 350;
% end;

Ns = size(G, 2) / 2; % two topography columns per each source of the grid
% perform projection of the coherence matrix away from the power only
Nch = size(G, 1);
% component
% vecC = CT(:); % Vectorize first
vecC = CT;
% project for each potential source
NsHR = size(G, 2) / 2; % two topography columns per each source of the grid
range3 = 1:3;
fprintf('Projecting cross-spectral data away from the power only component\n');

A = zeros(Nch ^ 2, Ns * 3);
for i = 1:NsHR
    range_i = i * 2 - 1:i * 2;
    gi = G(:,range_i);
    v = gi(:,1) * gi(:,1)';
    V(:,1) = v(:);
    v = gi(:,2) * gi(:,2)';
    V(:,2) = v(:);
    v = gi(:,1) * gi(:,2)';
    V(:,3) = v(:);
    v = gi(:,2) * gi(:,1)';
    V(:,4) = v(:);
    [ans, sg0, vg0] = svd(V' * V, 'econ');
    vg0  = vg0(:,1:3) * diag((1. / sqrt(diag(sg0(1:3, 1:3)))));
    ug = (V * vg0);
    A(:, range3) = ug;
    range3 = range3 + 3;
end;
AA = A * A';
[Upwr, e] = eigs(AA, Rnk);
vecCp  = vecC - Upwr * Upwr' * vecC;
clear AA A;
Cp = vecCp;
if(~isempty(CT2))
    Rnk12 = 20;
    CT2  = CT2-Upwr * (Upwr' * CT2);
    [u2,s2,v2] = svd(CT2,'econ');
    Cp = Cp-u2(:,1:Rnk12)*(u2(:,1:Rnk12)'*Cp); 
end;
% put back in shape
% Cp = reshape(vecCp, size(C,1), size(C,2));
% initialize source space coupling indicator vector

