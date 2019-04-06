function [] = get_Q24_final_MSE

% This function obtains the MSE of model parameters under a 4*2 Q-matrix,
% corresponding to Figure 2 in Section 3.1

%%%%%%%%%%%%%%%%% Q24: random generation and generic identifiability %%%%%%%%%%%%%%%%%%%
Q = [1 0; 0 1; 1 0; 0 1];
[J, K] = size(Q);

% I_use = get_I(2, J);
% I_full = get_I(J,J);
% I_full1 = [zeros(1,J); I_full];
% 
% A = binary(0:(2^K-1), K);

% randomly generate Nset sets of parameters
Nset = 100;
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


[p_mse, c_mse, g_mse] = get_Q24_single_MSE(N, p_true, c_true, g_true);

mean(p_mse)
mean(c_mse)
mean(g_mse)

filename = strcat('New_Q24_MSE_N', num2str(N), '.mat');
save(filename)


end