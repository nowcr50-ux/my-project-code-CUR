function [C,U,R] = FastAdaCUR(A,t,e,b,p)

[~,n,~]= size(A);
q = length(t);

C = cell(1,q);
U = cell(1,q);
R = cell(1,q);

A1 = A(:,:,t(1));

[I,J] = Rand_Pivot_RankEst(A1,e/sqrt(n));
r = length(I);
I0 = OS(A1,I,J,p+b);
J0 = OS(A1',J,I,b);
%length(J0)
%I = [I I0];
J = [J J0];
I = unique([I I0], 'stable');
%J = unique([J J0], 'stable');



C{1}= A1(:,J(1:r));
R{1} =A1(I(1:r+p),:);
U{1} = A1(I(1:r+p),J(1:r));


for j =2:q
    Aj = A(:,:,t(j));
    Ue = Aj(I,J);
    %addition
     % if j==2
     %     p=10;
     % end
    %endaddition

     %[~,~,I1]=sRRQR(Ue',2,"tol",0.05*e/sqrt(n));
     %[~,Re,J1] =sRRQR(Ue,2,"tol",0.05*e/sqrt(n));
     [~,~,I1]=sRRQR(Ue',2,"rank",r+b+p);
     %I1= CPQR(Ue');
     [~,Re,J1] =sRRQR(Ue,2,"rank",r+b);

    d = abs(diag(Re));
    r0 = sum(d >(e/sqrt(n)) * abs(Re(1,1)));
    

    if r0<=r
        
       I = I(I1);
       J = J(J1);

       I = I(1:r0+b+p);
       %J = J(1:r0+b); original
       J = J(min(1:r0+b,length(J)));

       % if length(I)<r0+b+p
       %     fprintf("Line 46 FastAdaCUR, I length")
       % end
       % if length(J)<r0+b
       %     fprintf("Line 49 FastAdaCUR, J length")
       % end


    else
        I2 = OS(Aj,I,J,r0-r);
        J2 = OS(Aj',J,I,r0-r);

        
        I = [I(I1) I2];
        J = [J(J1) J2];

        if length(I)<r0+b+p
           l = abs(length(I)-r0+b+p);
           I2 = OS(Aj,I,J,r0-r+l);
           
           I = [I(I1) I2]; 
            
        end
        if length(J)< r0+b
            ll = abs(length(J)-r0+b);
            J2 = OS(Aj',J,I,r0-r+ll);
            J = [J(J1) J2];

            fprintf('Line 62 FastAda,|J|')
        end
        
    end
    r = r0;
    %addition
    if j==2 || j==3 
    [I,J] = Rand_Pivot_RankEst(Aj,e/sqrt(n));

    I0 = OS(Aj,I,J,p+b);
    J0 = OS(Aj',J,I,b);
    r = length(I);

    J = [J J0];
    I = unique([I I0], 'stable');
    C{j} = Aj(:,J(1:r));
    R{j} = Aj(I(1:r+p),:);
    U{j} = Aj(I(1:r+p),J(1:r));
    r=r0;
    end
    %endaddition (comment block)

    C{j} = Aj(:,J(1:r));
    R{j} = Aj(I(1:r+p),:);
    U{j} = Aj(I(1:r+p),J(1:r));
end

end



