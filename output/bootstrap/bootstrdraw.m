
figure;
hctx  = trisurf(Ctx.Faces,Ctx.Vertices(:,1),Ctx.Vertices(:,2),Ctx.Vertices(:,3),'FaceColor',[0.1,0.51,1], 'EdgeColor','none','FaceAlpha', 0.1);
hold on;

camlight left; lighting phong
camlight right; 
hold on;
 % plot3(XYZGen(:,1), XYZGen(:,2),XYZGen(:,3),'r');
cols = ['g','r','m','y','k','b'];

occ = zeros(size(A));
for i = 1:length(A)
	occ(i) = sum(A == A(i));
end
B = A(occ>0);
Result = [(mod(B,Nsites))', ((B - mod(B,Nsites)) / Nsites+ 1)'];
[Npairs, dummy] = size(Result);
adjMat = zeros(Npairs,Npairs);
Dpair = 0.013; % Pairwise clustering distance threshold
for p1 = 1:Npairs
    for p2 = p1:Npairs
        if norm(R(Result(p1,1),:) - R(Result(p2,1),:)) < Dpair && norm(R(Result(p1,2),:) - R(Result(p2,2),:)) < Dpair
            adjMat(p1,p2) = 1;
            adjMat(p2,p1) = 1;
        end
    end
end

N = size(Result,1);
i = 1;
Pairs = Result;
while(~isempty(Pairs))
    if length(nonzeros(bfs(adjMat,i) > 0)) > 10
        % Pairs(i,:)
        clust = Pairs(bfs(adjMat,i) > -1,:);
        drawset(clust, R);
        
        Pairs = Pairs(bfs(adjMat,i) == -1,:);
        adjMat = adjMat(bfs(adjMat,i) == -1, bfs(adjMat,i) == -1);
    	% L1 = line( R(Pairs(i,:), 1), R(Pairs(i,:),2), R(Pairs(i,:),3) );
    	% plot3( R(Pairs(i,1), 1), R(Pairs(i,1),2), R(Pairs(i,1),3), 'c.', 'Markersize', 40 );
     %    plot3( R(Pairs(i,2), 1), R(Pairs(i,2),2), R(Pairs(i,2),3), 'c.', 'Markersize', 40 );
    	% set(L1, 'Color', 'c', 'linewidth',2);
        % i = 1;
    else
        Pairs = Pairs(2:end,:);
        adjMat = adjMat(2:end,2:end);
        % i = i + 1;
    end
end
hold on;


% i = DipInd(1);
% j = DipInd(2);
% range_i = i*2-1:i*2;
% ai = G2dU(:,range_i)';
% range_j = j*2-1:j*2;
% aj = G2dU(:,range_j)';
% cs = ai*real(Cp{1})*aj';
% [u s v] = svd(cs);
% a = u(:,1)'*ai;
% b = v(:,1)'*aj;
% abRe = a(:)*b(:)'+ b(:)*a(:)';
% abRe = abRe(:);
% cs = ai*imag(Cp{1})*aj';
% [u s v] = svd(cs);
% a = u(:,1)'*ai;
% b = v(:,1)'*aj;
% abIm = a(:)*b(:)'- b(:)*a(:)';
% abIm = abIm(:);
