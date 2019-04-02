function [p_true, c_true, g_true] = get_rand_true(Nset)

Q = [1 0; 0 1; 1 0; 0 1];
[J, K] = size(Q);

p_true = zeros(2^K, Nset);
c_true = zeros(J, Nset);
g_true = zeros(J, Nset);

rng('default');
for rr = 1:Nset
    % pn = 0.1 * rand(4,1) + 0.1;
    pn = gamrnd(ones(1,2^K), ones(1,2^K));
    p_true(:,rr) = pn / sum(pn);
    c_true(:,rr) = 0.9 - 0.3 * rand(J, 1);
    g_true(:,rr) = 0.1 + 0.3 * rand(J, 1);
end
    



end