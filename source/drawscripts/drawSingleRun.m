%% function drawSingleRun
%  drawSingleRun(A, R, col)
% A - set of numbers of pairs which want to draw 
% R - matrix of coordinates of grid nodes locations
% col - color specification
function drawSingleRun(A, R, col, Ctx)
 % col = randint(1,3);
 figure;
hctx  = trisurf(Ctx.Faces,Ctx.Vertices(:,1),Ctx.Vertices(:,2),Ctx.Vertices(:,3),'FaceColor',[0.1,0.51,1], 'EdgeColor','none','FaceAlpha', 0.1);
hold on;
Nsites = size(R,1);
DipInd = [(mod(A,Nsites))', ((A - mod(A,Nsites)) / Nsites+ 1)'];
for i = 1:length(DipInd)
		L1 = line( R(DipInd(i,:), 1), R(DipInd(i,:),2), R(DipInd(i,:),3) );
		plot3( R(DipInd(i,1), 1), R(DipInd(i,1),2), R(DipInd(i,1),3),'.', 'Color', col, 'Markersize', 40 );
        plot3( R(DipInd(i,2), 1), R(DipInd(i,2),2), R(DipInd(i,2),3),'.', 'Color', col, 'Markersize', 40 );
     	set(L1, 'Color', col, 'linewidth',2);
end