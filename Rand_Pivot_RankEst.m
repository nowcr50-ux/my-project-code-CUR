function [I,J] = Rand_Pivot_RankEst(A, e)
s0 = 10;
s_step = 5;


[r_hat, X] = rand_rank_est(A, e, s0, s_step);



[I,J] = Rand_Pivot(A, X, r_hat);
end