
figure;
hctx  = trisurf(Ctx.Faces,Ctx.Vertices(:,1),Ctx.Vertices(:,2),Ctx.Vertices(:,3),'FaceColor',[0.5,0.51, 0.5], 'EdgeColor','none','FaceAlpha', 0.1);
hold on;

camlight left; lighting phong
camlight right; 
hold on;
 % plot3(XYZGen(:,1), XYZGen(:,2),XYZGen(:,3),'r');
cols = ['g','r','m','y','k','b','c', 'w'];

occ = zeros(size(A));
for i = 1:length(A)
    occ(i) = sum(A == A(i));
end
B = A(occ > 0);
[Npairs, dummy] = size(DipInd);
adjMat = zeros(Npairs, Npairs);
Dpair = 0.0112; % Pairwise clustering distance threshold
clustSize = 20; % Clulstersize threshold 
for p1 = 1:Npairs
    for p2 = p1:Npairs
        if norm(R(DipInd(p1, 1),:) - R(DipInd(p2,1),:)) < Dpair ...
&& norm(R(DipInd(p1,2),:) - R(DipInd(p2,2),:)) < Dpair
            adjMat(p1,p2) = 1;
            adjMat(p2,p1) = 1;
        end
    end
end

N = size(DipInd, 1);
i = 1; % Number of column in adjacence martrix
clustNum = 0;
clusters = cell(1, 20);
restPairs = DipInd;
while(~isempty(restPairs))
    if length(nonzeros(bfs(adjMat,i) > 0)) > clustSize;
        % restPairs(i,:)
        clust = restPairs(bfs(adjMat,i) > -1,:);
        clustNum = clustNum + 1;
        drawset(clust, R, cols(clustNum));
        clusters{clustNum} = clust;
        restPairs = restPairs(bfs(adjMat, i) == -1,:);
        adjMat = adjMat(bfs(adjMat, i) == -1, bfs(adjMat, i) == -1);
        % L1 = line( R(restPairs(i,:), 1), R(restPairs(i,:),2), R(restPairs(i,:),3) );
        % plot3( R(restPairs(i,1), 1), R(restPairs(i,1),2), R(restPairs(i,1),3), 'c.', 'Markersize', 40 );
        % plot3( R(restPairs(i,2), 1), R(restPairs(i,2),2), R(restPairs(i,2),3), 'c.', 'Markersize', 40 );
        % set(L1, 'Color', 'c', 'linewidth',2);
        i = 1;
    else
        restPairs = restPairs(2:end,:);
        adjMat = adjMat(2:end,2:end);
        i = i + 1;
    end
end
hold on;

return;
[Nch, Nsrc] = size(G2dLRU);
for iClust = 1:clustNum
figure;
currClust = clusters{iClust};
plotTs = [];
for k = 1:numfiles
    Cp  = reshape(mean(M_res{k},2), Nch, Nch);
    for num = 1:length(mydata{k}.A)
        currPair = [(mod(mydata{k}.A(num),Nsites))', ((mydata{k}.A(num) - mod(mydata{k}.A(num), Nsites)) / Nsites+ 1)'];
        if( ~isempty(intersect(currPair, currClust, 'rows')) )
            i = currPair(1);
            j = currPair(2);
            range_i = i * 2-1:i * 2;
            ai = G2dLRU(:,range_i)';
            range_j = j * 2-1:j * 2;
            aj = G2dLRU(:,range_j)';
            cs = ai * real(Cp) * aj';
            [u s v] = svd(cs);
            a = u(:,1)' * ai;
            b = v(:,1)' * aj;
            abRe = a(:) * b(:)'+ b(:) * a(:)';
            abRe = abRe(:);
            cs = ai * imag(Cp) * aj';
            [u s v] = svd(cs);
            a = u(:,1)' * ai;
            b = v(:,1)' * aj;
            abIm = a(:) * b(:)'- b(:) * a(:)';
            abIm = abIm(:);

            CpTime  = M_res{k}; % projected away from power time-resolved cross-spectrum
             
            TsRe = abRe(:)' * real(CpTime);
            TsIm = abIm(:)' * imag(CpTime);

            Ts = sqrt(TsRe .^ 2 + TsIm .^ 2);
            plotTs = [plotTs; Ts];
            % plot(1:351, Ts);
            % hold on;
        end
    end

end
plot(1:351, mean(plotTs, 1), 'linewidth', 2, 'Color', cols(iClust));
% errorbar(1:351, mean(plotTs,1), std(plotTs,1));
end
