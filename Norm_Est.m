function E = Norm_Est(A,I,J,s)
[m,~] = size(A);
Gam = randn(s,m);

X = Gam*A;

C = A(:,J);
U = A(I,J);
R = A(I,:);

X1 = Gam * (C *pinv(U)*R );


E = norm (X-X1,"fro");
end




