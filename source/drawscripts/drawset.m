%% function drawset(A, R, col)
% Draws set A using grid locations R and color col
function drawset(A, R, col)
 % col = randint(1,3);
    linewidth = 1;
    Markersize = 40;
    for i = 1:length(A)
            L1 = line( R(A(i,:), 1), R(A(i,:),2), R(A(i,:),3) );
            plot3( R(A(i,1), 1), R(A(i,1),2), R(A(i,1),3),'.', 'Color', col, 'Markersize', Markersize );
            plot3( R(A(i,2), 1), R(A(i,2),2), R(A(i,2),3),'.', 'Color', col, 'Markersize', Markersize );
            set(L1, 'Color', col, 'linewidth', linewidth);
    end
