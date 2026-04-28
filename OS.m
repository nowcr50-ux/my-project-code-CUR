function I0 = OS(A,I,J,p)
[m,~]=size(A);
r = length(I);

       
[Qc,~]=qr(A(:,J));
[~,~,V]=svd(Qc(I,:));
      
p = min(p, r); 
if p == 0
   I0 = [];
   return;
end
Vp = V(:,r-p+1:r);

K = setdiff(1:m, I, 'stable');

if isempty(K)
    I0 = [];
    return;
end

QC = Qc(K,:)*Vp;
i = CPQR(QC');
%I0 = i(1:p);
I0 = K(i(1:p));

end

