%% columnG: returns column of G; p - number of column, G_small - generating model matrix in physical space.
function [column] = columnG(p, G_small)
	[Nsen, Nsrc] = size(G_small); % N equals number of sources in physical space, Sen - num. of sensors in phys. space
							  % (opposed to cross-correlation space, which comprises N^2 sources)
	i = mod(p, Nsrc);
	if i == 0
		i = Nsrc;
	end
	j = (p - i) / Nsrc + 1;
	column(Nsen^2) = 0;
	matColumn = G_small(:,i) * G_small(:,j)';
	column = matColumn(:);
