%%% Construct examples with K=3 to verify the necessary condition for 
%%% generic identifiability of GDINA with unknown Q %%%

% %%%%%% one true Q %%%%%%%
% Q = [1 0 0; 1 0 0; [zeros(2,1), eye(2)]; [zeros(2,1), eye(2)]; 0 1 1];
% 
% Q_bar = [1 1 1; 1 1 1; [zeros(2,1), eye(2)]; [zeros(2,1), eye(2)]; 0 1 1];


%%%%%% another true Q %%%%%%%
Q = [1 1 0; 1 0 1; [zeros(2,1), eye(2)]; [zeros(2,1), eye(2)]; 0 1 1];

Q_bar = [1 1 1; 1 1 1; [zeros(2,1), eye(2)]; [zeros(2,1), eye(2)]; 0 1 1];


[J, K] = size(Q);

I_use = get_I(2, J);
I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

nu_true = ones(2^K, 1)/2^K;

p1 = 0.2; p2 = 0.8;
delta_true = [p1*ones(J,1), zeros(J,2^K-1)];
theta_true = zeros(J, 2^K);


S = sum(Q, 2);
expand_S = cell(K, 1);
for j = 1:J
    expand_S{j} = [zeros(1, S(j)); get_I(S(j), S(j))];
end


attr_combo = get_I(K, K); % (2^K-1) * K
expand_Q = prod(bsxfun(@power, reshape(Q, [J 1 K]), reshape(attr_combo, [1 2^K-1 K])), 3);
expand_Q = [ones(J,1), expand_Q];
for j = 1:J
    theta_temp = zeros(2^S(j), 1);
    delta_temp = [p1; zeros(2^S(j)-1, 1)];
    delta_temp(2:2^S(j)) = (p2-p1)/(2^S(j)-1);
    
    step_prob = (p2-p1)/S(j);
    count = 1;
    for jj = 1:S(j)
        theta_temp((count+1):(count+nchoosek(S(j), jj))) = p1 + jj * step_prob;
        count = count + nchoosek(S(j), jj);
    end
    
    ok = (expand_Q(j,:)==1);
    delta_true(j, ok) = delta_temp(1:2^S(j));

    Aj = [zeros(1, S(j)); get_I(S(j), S(j))];
    % expand Aj to get Mj
    Mj = prod(bsxfun(@power, reshape(Aj, [2^S(j) 1 S(j)]), reshape(Aj, [1 2^S(j) S(j)])), 3);
    theta_true(j, ok) = (Mj * delta_true(j, ok)')';
end

attr_combo = get_I(K, K);
A = binary(0:(2^K-1), K);
expand_A = [ones(size(A,1), 1), prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];

prob_posi_true = delta_true * expand_A';
T_use_true = get_T_GDINA(Q, I_use, prob_posi_true);
PR_cond_true = get_rp_GDINA(Q, I_full, prob_posi_true);
PR_marg_true = PR_cond_true * nu_true;

Theta_true = T_use_true(1:J, :);



% Theta_bar(1, 1:4) = [0.18, 0.19, 0.21, 0.22];
% Theta_bar(2, 1:4) = [0.23, 0.24, 0.25, 0.26];


%%%%% begin trial %%%%%
Theta_bar = Theta_true; nu_bar = zeros(1, 2^K);

% rng(1); rng(30); rng(100); rng(1300); rng(1400); rng(1700)
rng(1700)
Theta_bar(1:2, 1:4) = Theta_true(1:2, 1:4) + (-0.1 + 0.2 * rand(2, 4));
Theta_bar(1, 5:8) = Theta_true(1, 1:4) + (Theta_true(1, 5:8) - Theta_true(1, 1:4)) .* (Theta_true(2,5:8) - Theta_bar(2,1:4)) .* (nu_true(5:8))' ./ ...
    ((Theta_true(2,1:4) - Theta_bar(2,1:4)) .* (nu_true(1:4))' + (Theta_true(2,5:8) - Theta_bar(2,1:4)) .* (nu_true(5:8))');

Theta_bar(2, 5:8) = Theta_true(2, 1:4) + (Theta_true(2, 5:8) - Theta_true(2, 1:4)) .* (Theta_true(1,5:8) - Theta_bar(1,1:4)) .* (nu_true(5:8))' ./ ...
    ((Theta_true(1,1:4) - Theta_bar(1,1:4)) .* (nu_true(1:4))' + (Theta_true(1,5:8) - Theta_bar(1,1:4)) .* (nu_true(5:8))');

nu_bar(5:8) = ( (Theta_true(2,1:4) - Theta_bar(2,1:4)) .* nu_true(1:4)' + (Theta_true(2,5:8) - Theta_bar(2,1:4)) .* nu_true(5:8)' )  ./ ...
    (Theta_bar(2,5:8) - Theta_bar(2,1:4));
nu_bar(1:4) = nu_true(1:4)' + nu_true(5:8)' - nu_bar(5:8);
nu_bar = nu_bar';

PR_cond_bar = get_rp_GDINA(Q, I_full, Theta_bar);
PR_marg_bar = PR_cond_bar * nu_bar;
max(abs(PR_marg_true - PR_marg_bar))

Theta_bar
nu_bar'
Theta_bar6 = Theta_bar; nu_bar6 = nu_bar;

% % check four equations
% nu_true(1:4)'+nu_true(5:8)' - nu_bar(1:4)'-nu_bar(5:8)'
% Theta_true(1,1:4) .* nu_true(1:4)' + Theta_true(1,5:8) .* nu_true(5:8)' - Theta_bar(1,1:4) .* nu_bar(1:4)' - Theta_bar(1,5:8) .* nu_bar(5:8)'
% Theta_true(2,1:4) .* nu_true(1:4)' + Theta_true(2,5:8) .* nu_true(5:8)' - Theta_bar(2,1:4) .* nu_bar(1:4)' - Theta_bar(2,5:8) .* nu_bar(5:8)'
% 
% Theta_true(1,1:4) .* Theta_true(2,1:4) .* nu_true(1:4)' + Theta_true(1,5:8) .* Theta_true(2,5:8) .* nu_true(5:8)' -...
%     Theta_bar(1,1:4) .* Theta_bar(2,1:4) .* nu_bar(1:4)' - Theta_bar(1,5:8) .* Theta_bar(2,5:8) .* nu_bar(5:8)'




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(1:2^J, PR_marg_true)
hold on
plot(1:2^J, PR_marg_bar)



figure
plot(1:16, [Theta_true(1,:), Theta_true(2,:)], '-ko', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:16, [Theta_bar1(1,:), Theta_bar1(2,:)], ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:16, [Theta_bar2(1,:), Theta_bar2(2,:)], ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:16, [Theta_bar3(1,:), Theta_bar3(2,:)], ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:16, [Theta_bar4(1,:), Theta_bar4(2,:)], ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:16, [Theta_bar5(1,:), Theta_bar5(2,:)], ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:16, [Theta_bar6(1,:), Theta_bar6(2,:)], ':+', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Item Parameters')
% lgd = legend('True', 'Theta^1', 'Theta^2', 'Theta^3', 'Location', 'northwest');
% lgd.FontSize = 12;
axis([1 16 0.1 1])
pbaspect([4 3 1]);
xticks(1:16)
xticklabels({'\theta_{1,000}','\theta_{1,001}','\theta_{1,010}','\theta_{1,011}',...
    '\theta_{1,100}','\theta_{1,101}', '\theta_{1,110}', '\theta_{1,111}', '\theta_{2,000}', '\theta_{2,001}',...
    '\theta_{2,010}', '\theta_{2,011}', '\theta_{2,100}', '\theta_{2,101}', '\theta_{2,110}', '\theta_{2,111}'})
xtickangle(45)
set(gca,'FontSize',12)
print('-r500', 'id_Q_Theta', '-dpng');



figure
plot(1:8, nu_true, '-ko', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:8, nu_bar1, ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:8, nu_bar2, ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:8, nu_bar3, ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:8, nu_bar4, ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:8, nu_bar5, ':+', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(1:8, nu_bar6, ':+', 'LineWidth', 2, 'MarkerSize', 8)
grid on
title('Proportion Parameters')
% lgd = legend('True', 'p^1', 'p^2', 'p^3', 'Location', 'northwest');
% lgd.FontSize = 12;
xticks(1:16)
xticklabels({'p_{000}','p_{001}','p_{010}','p_{011}','p_{100}',...
    'p_{101}','p_{110}', 'p_{111}'})
%xtickangle(45)
pbaspect([4 3 1]);
set(gca,'fontsize',14)
print('-r500', 'id_Q_p', '-dpng');
