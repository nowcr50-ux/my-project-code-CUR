function [Is,Js] = Pivot_Residual(A,XR,I,J)
[s,~]= size(XR);
Js = CPQR(XR);
Js = Js(~ismember(Js, J)); 
Js = Js(1:s); 

C = A(:,J);
U = A(I,J);
Rjs = A(I,Js);

Cr = C *(U\Rjs);
CR = A(:,Js) - Cr;

Is = CPQR(CR');
Is = Is(~ismember(Is, I)); 
Is = Is(1:s); 
end
