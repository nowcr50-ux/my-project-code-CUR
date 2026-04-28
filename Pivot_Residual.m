function [Is,Js] = Pivot_Residual(A,XR,I,J)
[s,~]= size(XR);
Js = CPQR(XR);
Js = Js(~ismember(Js, J)); %added
Js = Js(1:s); %original

C = A(:,J);
U = A(I,J);
Rjs = A(I,Js);

Cr = C *(U\Rjs);
CR = A(:,Js) - Cr;

Is = CPQR(CR');
Is = Is(~ismember(Is, I)); %added
Is = Is(1:s); %original
end