function [I,J] = Rand_Pivot(A,X,r)
Jj = CPQR(X);
J = Jj(1:r);
L = A(:,J)';
Ii = CPQR(L);
I= Ii(1:r);

end