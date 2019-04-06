% This file corresponds to Simulation Study VII in the Supplementary
% Material for the case K=5. This file generates Figure 9 in the supplement

% K = 5, J = 20
Q_18times4 = [eye(4); eye(4); eye(4); 1 1 0 0; 1 0 1 0;...
    1 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1];

Q = [1 1 0 0 0; 1 0 1 0 0; [zeros(18,1), Q_18times4]];

Q_bar = Q;
Q_bar(1:2, :) = ones(2, 5);


% % 
[J, K] = size(Q);

A = binary(0:(2^K-1), K);

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


% get delta matrix from theta matrix
attr_combo = get_I(K, K); % (2^K-1) * K
expand_Q = [ones(J,1), prod(bsxfun(@power, reshape(Q, [J 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];

ok_ex_all = (expand_Q==1);

MM = cell(J, 1);
inv_MM = cell(J, 1);
for j = 1:J
    MM{j} = prod(bsxfun(@power, reshape(expand_S{j}, ...
        [2^S(j) 1 S(j)]), reshape(expand_S{j}, [1 2^S(j) S(j)])), 3);
    inv_MM{j} = inv(MM{j});
end

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


expand_A = [ones(size(A,1), 1), prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];

% columns of prob_posi are arranged in order as rows of A (binary order)
prob_posi_true = delta_true * expand_A';
T_use_true = get_T_GDINA(Q, I_use, prob_posi_true);
PR_cond_true = get_rp_GDINA(Q, I_full, prob_posi_true);
PR_marg_true = PR_cond_true * nu_true;

% columns of Theta_true are arranged in order as rows of A (binary order)
Theta_true = T_use_true(1:J, :);



%%%%% begin trial: generate parameters under Q_bar %%%%%
% first expand Q_bar
expand_Qbar = [ones(J,1), prod(bsxfun(@power, reshape(Q_bar, [J 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];

num_try=1000;
is_mono_Qbar = zeros(num_try,1); 
is_in_range = zeros(num_try,1);
is_rough_id = zeros(num_try,1);

delta_bar_arr = zeros(J, 2^K, num_try); 
theta_bar_arr = zeros(J, 2^K, num_try);
nu_bar_arr = zeros(2^K, 1, num_try);
for ii = 1:num_try
    % rng(ii)
    [Theta_bar, theta_bar_combo, delta_bar_combo, nu_bar, is_monotone] = ...
        get_gdina_theta_bar_THM3(Theta_true, nu_true, Q_bar, expand_Qbar, ii);
    
    is_mono_Qbar(ii) = is_monotone;
    delta_bar_arr(:,:,ii) = delta_bar_combo;
    theta_bar_arr(:,:,ii) = theta_bar_combo;
    nu_bar_arr(:,:,ii) = nu_bar;
    
    is_in_range(ii) = all(all(Theta_bar<=1 & Theta_bar>=0));
    
    is_rough_id(ii) = all(max(Theta_bar(:,end) > Theta_bar(:,1:(end-1)),[],2));
    
    fprintf('Run %d completed, parameter within range %d, status of rough monotonicity: %d,  \n',...
        ii, is_in_range(ii), is_rough_id(ii));
end

% find_rng0 = find(is_mono_Qbar==1 & is_in_range==1)

find_rng0 = find(is_rough_id==1 & is_in_range==1)

find_rng = find_rng0(1:70);



% Plot the parameter values to contrast their difference
Theta_bar_mono = zeros(J, 2^K, length(find_rng));
nu_bar_mono = zeros(2^K, 1, length(find_rng));

% load('nolid_gdina_attr2_K5_J20_70alter.mat')

% plot the item parameters
figure
for ind = 1:length(find_rng)
    [Theta_bar, theta_bar_combo, delta_bar_combo, nu_bar, is_monotone] = ...
        get_gdina_theta_bar_THM3(Theta_true, nu_true, Q_bar, expand_Qbar, find_rng(ind));
    
    Theta_bar_mono(:,:,ind) = Theta_bar;
    nu_bar_mono(:,:,ind) = nu_bar;
    plot(1:(2*2^K), [Theta_bar(1,:), Theta_bar(2,:)]-[Theta_true(1,:), Theta_true(2,:)],...
        ':+', 'LineWidth', 1, 'MarkerSize', 5)
    hold on
end
plot(1:(2*2^K), [Theta_true(1,:), Theta_true(2,:)]-[Theta_true(1,:), Theta_true(2,:)],...
    '-ko', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 3)
% title('Item Parameters')
% axis([1 (2*2^K) 0.1 1])
xlim([1 2*2^K])
pbaspect([4 3 1]);
set(gca,'FontSize',13)
xlabel('indices of item parameters')
ylabel('difference between alternative and true parameters')
grid on
print('-r300', 'nolid_gdina_attr2_K5_itemdif', '-dpng');


% plot the proportion parameters
figure
for ind = 1:length(find_rng)
    nu_bar = nu_bar_mono(:,:,ind);
    plot(1:(2^K), nu_bar'-nu_true', ':+', 'LineWidth',  1, 'MarkerSize', 5)
    hold on
end
plot(1:(2^K), nu_true'-nu_true', '-ko', 'MarkerFaceColor', 'k', 'LineWidth',  1, 'MarkerSize', 5)
% title('Item Parameters')
% axis([1 2^K 0 0.5])
xlim([1 2^K])
pbaspect([4 3 1]);
set(gca,'FontSize',14)
xlabel('indices of proportion parameters')
ylabel('difference between alternative and true parameters')
grid on
print('-r300', 'nolid_gdina_attr2_K5_prop', '-dpng');

% save('nolid_gdina_K3_J20_attr2_70alter.mat')

% For those alter. sets satisfying monotone constraint, check difference of
% marginald distributions of R
max_diff = zeros(length(find_rng), 1);
for ind = 1:length(find_rng)
    
    PR_cond_bar = get_rp_GDINA(Q_bar, I_full, Theta_bar_mono(:,:,ind));
    PR_marg_bar = PR_cond_bar * nu_bar_mono(:,:,ind);
    max_diff(ind) = max(abs(PR_marg_true - PR_marg_bar));
    fprintf('Run %d completed\n', ind);
end
max(max_diff)


