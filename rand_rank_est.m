function [r_hat, X] = rand_rank_est(A, eps_tol, s0, s_step)
    [m, n] = size(A);
    s = s0;

    while true

 
        Gamma1 = randn(s, m)/sqrt(s);
        Gamma2 = randn(n, 2*s)/sqrt(2*s);

      
        X = Gamma1 * A;     
        B = X * Gamma2;     

      
        sv = svd(B);

        
        normA_est = sv(1);
        r_hat = sum(sv > eps_tol * normA_est);

        r_hat = min(r_hat, min(size(A))); %added this line

        if sv(end) < eps_tol * normA_est
             break
        else
         s = s + s_step;
        end
    end
end