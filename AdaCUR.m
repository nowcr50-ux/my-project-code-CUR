function [C,U,R] = AdaCUR(A,t,e,s,p)

[m,n,~]= size(A);
q = length(t);

minorcount =0;
scratchcount =1;



C = cell(1,q);
U = cell(1,q);
R = cell(1,q);

A1 = A(:,:,t(1)); 

[I,J] = Rand_Pivot_RankEst(A1,0.5*e/sqrt(n));


I0 = OS(A1,I,J,p);
%fprintf("size of I0 is =%d\n",length(I0));

In = [I I0]; 


C{1} = A1(:,J);
R{1} = A1(In,:);
U{1} = A1(In,J);



for j =2:q
    Aj = A(:,:,t(j));     
    gam = randn(s,m);
    X = gam * Aj;
  
    Cj = Aj(:,J);
    Uj = Aj(In,J);
    Rj = Aj(In,:);

    L = Cj * (Uj\Rj);
    Lj = gam*L;
    Es = X - Lj;

    E = norm(Es,"fro")/norm(X,"fro");
  
    

    if E> e
        [I1, J1] = Pivot_Residual(Aj, Es, I, J); 
        I=[I I1];
        J =[J J1];   
        In = [I I0];
     
        Ar = Aj(In,J);
        
       
        [~,~,I2]= sRRQR(Ar',2,"tol",0.5*e/sqrt(n));
        [~,Re,J2] = sRRQR(Ar,2,"tol",0.5*e/sqrt(n));


        d = abs(diag(Re));
        r = sum(d > (e/sqrt(n)) * abs(Re(1,1)));
        I = In(I2(1:r));
        %In = [I I0];
        % if r+p> length(I2)
        %     fprintf("Line 64 Ada \n");
        % end
        I0= In(I2(r+1:min(r+p,length(I2))));
        In = [I I0]; 
        J = J(J2(1:r));

        Cj = Aj(:,J);
        Uj = Aj(In,J);
        Rj = Aj(In,:);

        minorcount = minorcount +1;


        L = Cj * (Uj\Rj);
        Lj = gam*L;
        Es = X -Lj;
        E = norm(Es,"fro")/norm(X,"fro");
        if E>e
            Aj= A(:,:,t(j)); 
            [I,J] = Rand_Pivot_RankEst(Aj,0.5*e/sqrt(n));
            I0 = OS(Aj,I,J,p);
            In = [I I0]; 
            %minorcount = minorcount-1;
            scratchcount = scratchcount+1;
        end

    end
    C{j} = Aj(:,J);
    R{j} = Aj(In,:);
    U{j}= Aj(In,J);
end
fprintf('Minor changes=%d and count from scratch=%d\n',minorcount,scratchcount)
end