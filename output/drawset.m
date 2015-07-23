function res = drawset(A, R)
 col = randint(1,3);
for i = 1:size(A,1)
		L1 = line( R(A(i,:), 1), R(A(i,:),2), R(A(i,:),3) );
		plot3( R(A(i,1), 1), R(A(i,1),2), R(A(i,1),3),'.', 'Color', col, 'Markersize', 40 );
        plot3( R(A(i,2), 1), R(A(i,2),2), R(A(i,2),3),'.', 'Color', col, 'Markersize', 40 );
    	set(L1, 'Color', col, 'linewidth',2);
end