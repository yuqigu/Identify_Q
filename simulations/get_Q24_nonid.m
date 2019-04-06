function [c, g, p, PR_marg_bar, max_diff] = get_Q24_nonid(ratio1_bar, ratio2_bar, g1_bar, g2_bar)

% This function calculates alternative sets of parameters under Q_{4*2} 
% corresponding to Figure 1 in Section 3.1
% When p_{00}p_{11} = p_{01}p_{10}, the parameters are not identifiable


Q = [1 0; 0 1; 1 0; 0 1];
[J, K] = size(Q);

I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

A = binary(0:(2^K-1), K);

p_true = 0.25 * ones(4, 1);
c_true = 0.8 * ones(J, 1);
g_true = 0.2 * ones(J, 1);


ratio1 = p_true(3)/p_true(1);
ratio2 = p_true(2)/p_true(1);

p = zeros(2^K, 1); g = zeros(J, 1); c = zeros(J, 1);


p(1) = 1/(1 + ratio1_bar + ratio2_bar + ratio1_bar*ratio2_bar);

p(2) = ratio2_bar * p(1);
p(3) = ratio1_bar * p(1);
p(4) = ratio2_bar * ratio1_bar * p(1);

g(1) = g1_bar;  g(2) = g2_bar;

c(3) = (g_true(3)*(g_true(1)-g(1)) + ratio1 * c_true(3)*(c_true(1)-g(1))) ./...
    ((g_true(1)-g(1)) + ratio1 * (c_true(1)-g(1)));

c(4) = (g_true(4)*(g_true(2)-g(2)) + ratio2 * c_true(4)*(c_true(2)-g(2))) ./...
    ((g_true(2)-g(2)) + ratio2 * (c_true(2)-g(2)));


c(1) = g(1) + ((g_true(1)-g(1))*(p_true(1)+p_true(2)) + (c_true(1)-g(1))*(p_true(3)+p_true(4))) / ...
    (p(1) * ratio1_bar * (1+ratio2_bar));


c(2) = g(2) + ((g_true(2)-g(2))*(p_true(1)+p_true(3)) + (c_true(2)-g(2))*(p_true(2)+p_true(4))) / ...
    (p(1) * ratio2_bar * (1+ratio1_bar));


g(3) = (c_true(3)*(p_true(3)+p_true(4)) + g_true(3)*(p_true(1)+p_true(2)) - c(3)*(p(3)+p(4))) / ...
    (p(1)+p(2));

g(4) = (c_true(4)*(p_true(2)+p_true(4)) + g_true(4)*(p_true(1)+p_true(3)) - c(4)*(p(2)+p(4))) / ...
    (p(1)+p(3));




PR_cond_true = get_rp(A, Q, c_true, g_true, I_full1);
PR_marg_true = PR_cond_true * p_true;

PR_cond_bar = get_rp(A, Q, c, g, I_full1);
PR_marg_bar = PR_cond_bar * p;


max_diff = max(abs(PR_marg_true - PR_marg_bar));



end