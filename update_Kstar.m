function [Kstar] = update_Kstar(Z)

Z = (Z + Z')/2;
[U, D] = eig(Z); 
D(find(D<0)) = 0;
Kstar = U*D*U';

end

