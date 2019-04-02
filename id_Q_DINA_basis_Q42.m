%%% example 3 %%%
% Q =      [1 0 0; 0 1 0; 1 1 1; 1 0 1; 1 1 0; 1 0 1; 1 1 0];
% 
% Q_bar1 = [1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 1 0 1; 1 1 0];
% 
% Q_bar2 = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 1 1 0; 1 0 1; 1 1 0];

Q =      [1 0 0; 0 1 0; 1 1 1; 1 0 0; 1 1 0; 1 1 1];

Q_bar1 = [1 0 0; 0 1 0; 0 1 1; 1 0 0; 1 1 0; 1 1 1];

Q_bar2 = [1 0 0; 0 1 0; 0 0 1; 1 0 0; 1 1 0; 1 1 1];


[J, K] = size(Q);

I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

A = binary(0:(2^K-1), K);


% ideal_resp = prod( bsxfun(@power, reshape(A, [1 size(A,1) K]), reshape(Q, [J 1 K])), 3);
% ideal_resp1 = prod( bsxfun(@power, reshape(A, [1 size(A,1) K]), reshape(Q_bar1, [J 1 K])), 3);
% ideal_resp2 = prod( bsxfun(@power, reshape(A, [1 size(A,1) K]), reshape(Q_bar2, [J 1 K])), 3);
ideal_resp = get_ideal_resp(Q, A);
ideal_resp1 = get_ideal_resp(Q_bar1, A);
ideal_resp2 = get_ideal_resp(Q_bar2, A);


rng(1)
pn = 0.1 * rand(8,1) + 0.1;
p_true = pn / sum(pn);
c_true = 0.9 - 0.2 * rand(J, 1);
g_true = 0.1 + 0.2 * rand(J, 1);

PR_cond_Q = get_rp(A, Q, c_true, g_true, I_full1);
PR_marg_Q = PR_cond_Q * p_true;



p_1 = p_true;
p_1(4) = 0;
p_1(3) = p_true(3) + p_true(4);

PR_cond_Q1 = get_rp(A, Q_bar1, c_true, g_true, I_full1);
PR_marg_Q1 = PR_cond_Q1 * p_1;


max(abs(PR_marg_Q - PR_marg_Q1))


p_2 = p_true;
p_2(2) = 0;  p_2(1) = p_true(1) + p_true(2);
p_2(4) = 0;  p_2(3) = p_true(3) + p_true(4);
p_2(6) = 0;  p_2(5) = p_true(5) + p_true(6);

PR_cond_Q2 = get_rp(A, Q_bar2, c_true, g_true, I_full1);
PR_marg_Q2 = PR_cond_Q2 * p_2;

max(abs(PR_marg_Q - PR_marg_Q2))



figure
plot(PR_marg_Q, 'o', 'MarkerSize', 6)
hold on
plot(PR_marg_Q1, 'x', 'MarkerSize', 6)
hold on
plot(PR_marg_Q2, '+', 'Color', 'k', 'MarkerSize', 6)
grid on
pbaspect([16 9 1]);
xlabel('response pattern indices')
ylabel('marignal probabilities')
set(gca,'fontsize',12)
lgd = legend('True Q', 'Q^1', 'Q^2',  'Location', 'northeast');
lgd.FontSize = 12;
print('-r500', 'id_Q_basis', '-dpng');



%%%%%%%%%%%%%%%%%% Q24: nonidentifiability construction %%%%%%%%%%%%%%%%%%%%%
Q = [1 0; 0 1; 1 0; 0 1];
[J, K] = size(Q);

I_use = get_I(2, J);
I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

A = binary(0:(2^K-1), K);

rng(1)
p_true = 0.25 * ones(4, 1);
c_true = 0.8 * ones(J, 1);
g_true = 0.2 * ones(J, 1);

ratio1 = p_true(3)/p_true(1);
ratio2 = p_true(2)/p_true(1);

p = zeros(2^K, 1); g = zeros(J, 1); c = zeros(J, 1);


ratio1_bar = 0.8;
ratio2_bar = 1.2;
p(1) = 1/(1 + ratio1_bar + ratio2_bar + ratio1_bar*ratio2_bar);

p(2) = ratio2_bar * p(1);
p(3) = ratio1_bar * p(1);
p(4) = ratio2_bar * ratio1_bar * p(1);

g(1) = 0.18;  g(2) = 0.22;

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


max(abs(PR_marg_true - PR_marg_bar))


%%%% function %%%%%
[c2, g2, p2, PR_marg_bar2, max_diff2] = get_Q24(1.3, 0.7, 0.23, 0.19)


figure
plot(PR_marg_true, 'o', 'Color', 'k', 'MarkerSize', 15)
hold on
plot(PR_marg_bar, 'x', 'Color', [0    0.4470    0.7410],  'MarkerSize', 15)
hold on
plot(PR_marg_bar, '+', 'Color', [0.8500    0.3250    0.0980],   'MarkerSize', 15)
xticks(1:16)
xticklabels({'0000','1000','0100','0010','0001','1100','1010', '1001', '0110', '0101', '0011', '1110', '1101', '1011', '0111', '1111'})
xtickangle(45)
pbaspect([4 3 1]);
grid on
xlabel('response patterns')
ylabel('marignal probabilities')
set(gca,'fontsize',12)
lgd = legend('true model', 'alternative model 1', 'alternative model 2', 'Location', 'southwest');
lgd.FontSize = 12;
print('-r500', 'id_Q24_marg', '-dpng');


figure
plot([p_true; g_true; 1-c_true],  '-ko', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot([p; g; 1-c], ':s', 'LineWidth', 2, 'MarkerFaceColor', [0    0.4470    0.7410], 'MarkerSize', 8)
hold on
plot([p2; g2; 1-c2], ':s', 'LineWidth', 2, 'MarkerFaceColor', [0.8500    0.3250    0.0980], 'MarkerSize', 8)
grid on
xlabel('parameter indices')
ylabel('parameter values')
set(gca,'fontsize',12)
lgd = legend('true model', 'alternative model 1', 'alternative model 2', 'Location', 'southwest');

xlim([1 12])
xticks(1:12)
xticklabels({'p_{00}', 'p_{01}', 'p_{10}', 'p_{11}', 'g_1', 'g_2', 'g_3', 'g_4',  's_1', 's_2', 's_3', 's_4'})
grid on
print('-r500', 'id_Q24_para', '-dpng');



%%%%%%%%%%%%%%%%% Q24: random generation and generic identifiability %%%%%%%%%%%%%%%%%%%
Q = [1 0; 0 1; 1 0; 0 1];
[J, K] = size(Q);

I_use = get_I(2, J);
I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

A = binary(0:(2^K-1), K);

rng(1)
pn = 0.1 * rand(4,1) + 0.1;
p_true = pn / sum(pn);
c_true = 0.9 - 0.2 * rand(J, 1);
g_true = 0.1 + 0.2 * rand(J, 1);

rng(1)
[p_mse500, c_mse500, g_mse500] = get_MSE(500,p_true,c_true,g_true);
[p_mse100, c_mse100, g_mse100] = get_MSE(100,p_true,c_true,g_true);
[p_mse200, c_mse200, g_mse200] = get_MSE(200,p_true,c_true,g_true);
[p_mse300, c_mse300, g_mse300] = get_MSE(300,p_true,c_true,g_true);
[p_mse400, c_mse400, g_mse400] = get_MSE(400,p_true,c_true,g_true);


para_mat = [p_mse100, c_mse100, g_mse100;...
p_mse200, c_mse200, g_mse200;
p_mse300, c_mse300, g_mse300; 
p_mse400, c_mse400, g_mse400; 
p_mse500, c_mse500, g_mse500]; 

Nvec = 100:100:500;



% new Nvec
Nvec = 50:50:500;
para_mat = [];
for i = 1:10
    [p_mse0, c_mse0, g_mse0] = get_MSE(Nvec(i),p_true,c_true,g_true);
    para_mat = [para_mat; p_mse0, c_mse0, g_mse0];
end



figure
plot(Nvec(1:7), para_mat(1:7,1), '-o', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'Color', 'k', 'MarkerSize', 8)
% lgd = legend('MSE of p', 'Location', 'northeast');
% lgd.FontSize = 12;
lgd = legend('MSE of p', 'Location', 'northeast');
xticks(Nvec);
xlabel('N');
%ylabel('MSE');
% axis([50 350 0.008 0.018]);
set(gca,'FontSize',16)
pbaspect([1 1 1]);
grid on
print('-r400', 'p_mse', '-dpng');


figure
plot(Nvec(1:7), para_mat(1:7,2), '-s', 'MarkerFaceColor', [    0    0.4470    0.7410], 'LineWidth', 2, 'MarkerSize', 8)
lgd = legend( 'MSE of s', 'Location', 'northeast');
xticks(Nvec);
xlabel('N');
%ylabel('MSE');
% axis([50 350 0.014 0.026]);

set(gca,'FontSize',16)
pbaspect([1 1 1]);
grid on
print('-r400', 's_mse', '-dpng');


figure
plot(Nvec(1:7), para_mat(1:7,3), '-^', 'MarkerFaceColor', [0.8500    0.3250    0.0980], 'Color', [0.8500    0.3250    0.0980], 'LineWidth', 2, 'MarkerSize', 8)

lgd = legend( 'MSE of g', 'Location', 'northeast');
xticks(Nvec);
xlabel('N');
%ylabel('MSE');
%axis([50 500 0.014 0.026]);

set(gca,'FontSize',16)
pbaspect([1 1 1]);
grid on
print('-r400', 'g_mse', '-dpng');



