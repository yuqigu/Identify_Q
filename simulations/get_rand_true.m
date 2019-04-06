function [p_true, c_true, g_true] = get_rand_true(Nset)
% This function randomly generates a total number of Nset sets of model
% parameters under DINA for the 4*2 Q-matrix, following the generating
% mechanism described in simulation scenario (b) in Section 3.1

Q = [1 0; 0 1; 1 0; 0 1];
[J, K] = size(Q);

p_true = zeros(2^K, Nset);
c_true = zeros(J, Nset);
g_true = zeros(J, Nset);

rng('default');
for rr = 1:Nset
    pn = gamrnd(3*ones(1,2^K), ones(1,2^K));
    p_true(:,rr) = pn / sum(pn);
    c_true(:,rr) = 0.9 - 0.2 * rand(J, 1);
    g_true(:,rr) = 0.1 + 0.2 * rand(J, 1);
end
    



end