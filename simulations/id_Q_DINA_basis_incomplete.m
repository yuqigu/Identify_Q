% This file correponds to simulation study IV in the supplementary material


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First example: K=3, J=20
Qs = [1 0 0; 1 1 0; 1 1 1];

Q =      [1 0 0; 0 1 0; 1 1 1; repmat(Qs, [5 1]); 1 1 1; 1 1 1];

% Q_bar1 = [1 0 0; 0 1 0; 0 1 1; repmat(Qs, [5 1]); 1 1 1; 1 1 1];
% 
% Q_bar2 = [1 0 0; 0 1 0; 0 0 1; repmat(Qs, [5 1]); 1 1 1; 1 1 1];

Qs1 = [1 0 0; 1 1 0; 0 1 1];
Qs2 = [1 0 0; 1 1 0; 0 0 1];
Q_bar1 = [1 0 0; 0 1 0; 1 1 1; repmat(Qs1, [5 1]); 1 1 1; 1 1 1];
Q_bar2 = [1 0 0; 0 1 0; 1 1 1; repmat(Qs2, [5 1]); 1 1 1; 1 1 1];


[J, K] = size(Q);

I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

A = binary(0:(2^K-1), K);

% generate true parameters under Q
rng(1)
pn = 0.1 * rand(2^K,1) + 0.1;
p_true = pn / sum(pn);
c_true = 0.9 - 0.2 * rand(J, 1);
g_true = 0.1 + 0.2 * rand(J, 1);

PR_cond_Q = get_rp(A, Q, c_true, g_true, I_full1);
PR_marg_Q = PR_cond_Q * p_true;

% obtain the ideal response matrix (Gamma-matrix) of the Q-matrices
ideal_resp = get_ideal_resp(Q, A);
ideal_resp1 = get_ideal_resp(Q_bar1, A);
ideal_resp2 = get_ideal_resp(Q_bar2, A);

% construct 1st set of alternative parameters under Q_bar1 here
p_1 = p_true;
p_1(4) = 0;
p_1(3) = p_true(3) + p_true(4);
%
PR_cond_Q1 = get_rp(A, Q_bar1, c_true, g_true, I_full1);
PR_marg_Q1 = PR_cond_Q1 * p_1;
% check if the true param. and the 1st alter. set lead to same distribution
max(abs(PR_marg_Q - PR_marg_Q1))


% construct 2nd set of alternative parameters under Q_bar2 here
p_2 = p_true;
p_2(2) = 0;  p_2(1) = p_true(1) + p_true(2);
p_2(4) = 0;  p_2(3) = p_true(3) + p_true(4);
p_2(6) = 0;  p_2(5) = p_true(5) + p_true(6);
%
PR_cond_Q2 = get_rp(A, Q_bar2, c_true, g_true, I_full1);
PR_marg_Q2 = PR_cond_Q2 * p_2;
% check if the true param. and the 1st alter. set lead to same distribution
max(abs(PR_marg_Q - PR_marg_Q2))


% load('id_Q_DINA_basis_K3_J20.mat')
% load('id_Q_DINA_basis_K5_J20.mat')
% plot the marginal probabilities of the response patterns
figure
plot(PR_marg_Q, 'o', 'MarkerSize', 11)
hold on
plot(PR_marg_Q1, 'x', 'Color', [0.8500    0.3250    0.0980], 'MarkerSize', 10)
hold on
plot(PR_marg_Q2, '+', 'Color', [0.4660  0.6740  0.1880], 'MarkerSize', 10)
grid on
pbaspect([16 9 1]);
xlim([0 2^J]);  
%ylim([0.02*10^(-3) 0.2*10^(-3)]);
xlabel('response pattern indices')
ylabel('marignal probabilities')
set(gca,'fontsize',12)
lgd = legend('true model with $Q_1$', 'alternative $Q^{\prime}_1$', ...
    'alternative $Q^{\prime\prime}_1$', 'Location', 'northeast');
set(lgd,'Interpreter','latex');
lgd.FontSize = 12;

print('-r200', 'id_Q_basis_K3_J20', '-dpng');


% plot the proportion parameters
p_len = length(p_true);

% figure
% plot(p_true,  '-o', 'LineWidth', 2, 'MarkerSize', 11)
% hold on
% plot(p_1, ':x', 'LineWidth', 2,  'MarkerSize', 10)
% hold on
% plot(p_2, ':+', 'Color', [0.4660  0.6740  0.1880], 'LineWidth', 2,  'MarkerSize', 10)
% grid on
% xlabel('indices of 8 proportion parameters')
% ylabel('parameter values')
% set(gca,'fontsize',12)
% lgd = legend('true model', 'alternative model 1', 'alternative model 2', 'Location', 'northeast');
% lgd.FontSize = 12;
% xlim([1 p_len])
% xticks(1:p_len)
% xticklabels({'p_{000}', 'p_{001}', 'p_{010}', 'p_{011}', 'p_{100}', 'p_{101}', 'p_{110}', 'p_{111}'})
% pbaspect([16 9 1]);
% print('-r400', 'id_Q_basis_K3_J20_prop', '-dpng')

% % plot parameter differences %
figure
plot(p_true-p_true,  '-o', 'LineWidth', 2, 'MarkerSize', 11)
hold on
plot(p_1-p_true, ':x', 'Color', [0.8500    0.3250    0.0980], 'LineWidth', 2,  'MarkerSize', 10)
hold on
plot(p_2-p_true, ':+', 'Color', [0.4660  0.6740  0.1880], 'LineWidth', 2,  'MarkerSize', 10)
grid on
xlabel('indices of 8 proportion parameters')
ylabel('difference between alternative and true parameters')
set(gca,'fontsize',12)
lgd = legend('true model with $Q_1$', 'alternative $Q^{\prime}_1$',...
    'alternative $Q^{\prime\prime}_1$', 'Location', 'northeast');
set(lgd,'Interpreter','latex');
lgd.FontSize = 12;
xlim([1 p_len])
xticks(1:p_len)
xticklabels({'{000}', '{001}', '{010}', '{011}', '{100}', '{101}', '{110}', '{111}'})
pbaspect([16 9 1]);
print('-r200', 'id_Q_basis_K3_J20_propdif', '-dpng')




save('id_Q_DINA_basis_K3_J20.mat')

load('id_Q_DINA_basis_K5_J20.mat')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second example: K=5, J=20
Qs = [1 0 0 0 0; 1 1 0 0 0; 1 1 1 0 0; 1 1 1 1 0; 1 1 1 1 1];
Q = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 1 1 1 1 1; ...
    repmat(Qs, [3 1])];


% Q_bar2 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 1 1 1; ...
%           repmat(Qs, [3 1])];      
%       
% Q_bar4 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; ...
%           repmat(Qs, [3 1])];

Qs2 = [1 0 0 0 0; 1 1 0 0 0; 1 1 1 0 0; 1 1 1 1 0; 0 0 1 1 1];
Qs4 = [1 0 0 0 0; 1 1 0 0 0; 1 1 1 0 0; 1 1 1 1 0; 0 0 0 0 1];
Q_bar2 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 1 1 1; Qs2; ...
          repmat(Qs, [2 1])];      
Q_bar4 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; Qs4; ...
          repmat(Qs, [2 1])];
      
      

[J, K] = size(Q);

I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

A = binary(0:(2^K-1), K);

% generate true parameters under Q
rng(1)
pn = 0.1 * rand(2^K,1) + 0.1;
p_true = pn / sum(pn);
c_true = 0.9 - 0.2 * rand(J, 1);
g_true = 0.1 + 0.2 * rand(J, 1);

PR_cond_Q = get_rp(A, Q, c_true, g_true, I_full1);
PR_marg_Q = PR_cond_Q * p_true;


% obtain the ideal response matrix (Gamma-matrix) of the Q-matrices
ideal_resp = get_ideal_resp(Q, A);
ideal_resp2 = get_ideal_resp(Q_bar2, A);
ideal_resp4 = get_ideal_resp(Q_bar4, A);

% ideal_resp1 = get_ideal_resp(Q_bar1, A);
% ideal_resp3 = get_ideal_resp(Q_bar3, A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at how the proportions corresponding to Q_bar2 should be set
ind_bar2_diff = find(ideal_resp2 ~= ideal_resp);
% get the column indices (for attr. patterns) where IR and IR_bar2 differs
ind_bar2_col = (ind_bar2_diff - mod(ind_bar2_diff, J)) / J + 1;
% ideal_resp(:,[7 8 15 16 23 24])
% ideal_resp2(:,[7 8 15 16 23 24])

% construct 1st set of alternative parameters under Q_bar2 here
p_2 = p_true;
p_2(ind_bar2_col) = 0;
p_2(ind_bar2_col - 1) = p_true(ind_bar2_col - 1) + p_true(ind_bar2_col);

PR_cond_Q2 = get_rp(A, Q_bar2, c_true, g_true, I_full1);
PR_marg_Q2 = PR_cond_Q2 * p_2;

% check if the true param. and the 1st alter. set lead to same distribution
max(abs(PR_marg_Q - PR_marg_Q2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct 2nd set of alternative parameters here
% look at how the proportions corresponding to Q_bar2 should be set
ind_bar4_diff = find(ideal_resp4 ~= ideal_resp);
% get the column indices (for attr. patterns) where IR and IR_bar2 differs
ind_bar4_col = (ind_bar4_diff - mod(ind_bar4_diff, J)) / J + 1;
% ideal_resp(:, [1 2   3 4   5 6   7 8   9 10   11 12   13 14   15 16   17 18 ...
%     19 20   21 22   23 24   25 26   27 28   29 30])
% ideal_resp4(:, [1 2   3 4   5 6   7 8   9 10   11 12   13 14   15 16   17 18 ...
%     19 20   21 22   23 24   25 26   27 28   29 30])
p_4 = p_true;
p_4(ind_bar4_col) = 0;
p_4(ind_bar4_col - 1) = p_true(ind_bar4_col - 1) + p_true(ind_bar4_col);

PR_cond_Q4 = get_rp(A, Q_bar4, c_true, g_true, I_full1);
PR_marg_Q4 = PR_cond_Q4 * p_4;

% check if the true param. and the 2nd alter. set lead to same distribution
max(abs(PR_marg_Q - PR_marg_Q4))

save('id_Q_DINA_basis_K5_J20.mat')

load('id_Q_DINA_basis_K5_J20.mat')

% plot the marginal probabilities of the response patterns
figure
plot(PR_marg_Q, 'o', 'MarkerSize', 11)
hold on
plot(PR_marg_Q2, 'x', 'MarkerSize', 10)
hold on
plot(PR_marg_Q4, '+', 'Color', [0.4660  0.6740  0.1880], 'MarkerSize', 10)
grid on
pbaspect([16 9 1]);
ylim([0.02*10^(-3) 0.2*10^(-3)]);
xlim([0 2^J])
xlabel('response pattern indices')
ylabel('marignal probabilities')
set(gca,'fontsize',12)
lgd = legend('true model with $Q_2$', 'alternative $Q^{\prime}_2$', ...
    'alternative $Q^{\prime\prime}_2$', 'Location', 'northeast');
set(lgd,'Interpreter','latex');
lgd.FontSize = 12;
print('-r200', 'id_Q_basis_K5_J20_zoom', '-dpng');


% plot the proportion parameters
p_len = length(p_true);

figure
plot(p_true-p_true,  '-o', 'LineWidth', 2, 'MarkerSize', 11)
hold on
plot(p_2-p_true, ':x', 'LineWidth', 1,  'MarkerSize', 10)
hold on
plot(p_4-p_true, ':+', 'Color', [0.4660  0.6740  0.1880], 'LineWidth', 1,  'MarkerSize', 10)
grid on
xlabel('indices of 32 proportion parameters')
ylabel('difference between alternative and true parameters')
set(gca,'fontsize',12)
%axP = get(gca,'Position');
lgd = legend('true model with $Q_2$', 'alternative $Q^{\prime}_2$', ...
    'alternative $Q^{\prime\prime}_2$', 'Location', 'Best');
set(lgd,'Interpreter','latex');
%set(gca, 'Position', axP);
xlim([1 p_len])
xticks(1:p_len)
xticklabels({'{00000}', '{00001}', '{00010}', '{00011}', ...
    '{00100}', '{00101}', '{00110}', '{00111}', ...
    '{01000}', '{01001}', '{01010}', '{01011}', ...
    '{01100}', '{01101}', '{01110}', '{01111}', ...
    '{10000}', '{10001}', '{10010}', '{10011}', ...
    '{10100}', '{10101}', '{10110}', '{10111}', ...
    '{11000}', '{11001}', '{11010}', '{11011}', ...
    '{11100}', '{11101}', '{11110}', '{11111}'})
xtickangle(45)
pbaspect([16 9 1]);
print('-r200', 'id_Q_basis_K5_J20_propdif', '-dpng');


% [p_true'; p_2'; p_4']

