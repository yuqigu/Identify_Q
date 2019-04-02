%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GDINA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the index of the true Q-matrix used to generate the data
id0 = 29;

Q_aa = load('Q_aa.mat');
Q_arr = cell2mat(struct2cell(Q_aa));

vec = (1:121)';

Q0 = Q_arr(:,:,id0);
Q = Q0;


[J, K] = size(Q);

I_use = get_I(2, J);
I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

nu_true = ones(2^K, 1)/2^K;


S = sum(Q, 2);
expand_S = cell(K, 1);
for j = 1:J
    expand_S{j} = [zeros(1, S(j)); get_I(S(j), S(j))];
end


% % March 16: replace the following with a function generate_GDINA_theta(Q)
% p1 = 0.2; p2 = 0.8;
% delta_true = [p1*ones(J,1), zeros(J,2^K-1)];
% theta_true = zeros(J, 2^K);
% attr_combo = get_I(K, K); % (2^K-1) * K
% expand_Q = prod(bsxfun(@power, reshape(Q, [J 1 K]), reshape(attr_combo, [1 2^K-1 K])), 3);
% expand_Q = [ones(J,1), expand_Q];
% for j = 1:J
%     theta_temp = zeros(2^S(j), 1);
%     delta_temp = [p1; zeros(2^S(j)-1, 1)];
%     delta_temp(2:2^S(j)) = (p2-p1)/(2^S(j)-1);
%     
%     step_prob = (p2-p1)/S(j);
%     count = 1;
%     for jj = 1:S(j)
%         theta_temp((count+1):(count+nchoosek(S(j), jj))) = p1 + jj * step_prob;
%         count = count + nchoosek(S(j), jj);
%     end
%     
%     ok = (expand_Q(j,:)==1);
%     delta_true(j, ok) = delta_temp(1:2^S(j));
% 
%     Aj = [zeros(1, S(j)); get_I(S(j), S(j))];
%     % expand Aj to get Mj
%     Mj = prod(bsxfun(@power, reshape(Aj, [2^S(j) 1 S(j)]), reshape(Aj, [1 2^S(j) S(j)])), 3);
%     theta_true(j, ok) = (Mj * delta_true(j, ok)')';
% end


% attr_combo = get_I(K, K);
% A = binary(0:(2^K-1), K);
% expand_A = [ones(size(A,1), 1), prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), ...
%     reshape(attr_combo, [1 2^K-1 K])), 3)];
% 
% prob_posi_true = delta_true * expand_A';
% T_use_true = get_T_GDINA(Q, I_use, prob_posi_true);
% PR_cond_true = get_rp_GDINA(Q, I_full, prob_posi_true);
% PR_marg_true = PR_cond_true * nu_true;
% 
% Theta_true = T_use_true(1:J, :);

rng(123)

[theta_true, delta_true] = generate_GDINA_theta(Q);


rng(0)
pn = gamrnd(5 * ones(1,2^K), ones(1,2^K));
p_true = pn / sum(pn);

% set initial values
p_ini = p_true;
theta_ini = theta_true;

N = 10000;
[X, XA] = generate_X_GDINA(N, nu_true, Q, delta_true);

% set the number of random initializations
num_ran = 5;

p_arr = p_true;
theta_arr = theta_true;

ll = zeros(size(Q_arr, 3), 1);
for r = 1:size(Q_arr, 3)
    Q = Q_arr(:,:,r);
    
    % obtain expand_S for this Q-matrix
    S = sum(Q, 2);
    expand_S = cell(K, 1);
    for j = 1:J
        expand_S{j} = [zeros(1, S(j)); get_I(S(j), S(j))];
    end
    
    % Revision: 5 random initializations
    llk_vec = zeros(num_ran, 1);
    p_vec = zeros(2^K, num_ran);

    % create a 3-dimensional array to store the item parameters
    theta_vec = zeros(J, 2^K, num_ran);

%     [p_hat0, c_hat0, g_hat0, loglike0] = get_cg(X, Q_arr(:,:,r),...
%             p_true, c_true, c_true);
%         
%     p_vec(:,1) = p_hat0; c_vec(:,1)=c_hat0; g_vec(:,1)=g_hat0;
%     llk_vec(1) = loglike0;
        
    for ii = 1:num_ran
        pn0 = gamrnd(5 * ones(1,2^K), ones(1,2^K));
        p_ini = pn0 / sum(pn0);

        theta_ini = generate_GDINA_theta(Q);

        [theta_hat, ~, p_hat, loglik] = get_GDINAprob(X, Q, ...
            theta_ini, p_ini, expand_S);

        p_vec(:,ii) = p_hat; 
        theta_vec(:,:,ii) = theta_hat;

        llk_vec(ii) = loglik;
    end
    [~, ind] = max(llk_vec);
    
    p_arr(:, :, r) = p_vec(:,ind);
    theta_arr(:, :, r) = theta_vec(:,:,ind);

    ll(r) = llk_vec(ind);
    
    fprintf('Just finished iteration %d,\t with loglik %1.4f,\n', ...
        r, llk_vec(ind));
end

% ll = ll(1:14);

idx = id0; val = ll(idx);
[val_max, idx_max] = max(ll);

figure
fig_all = plot(setdiff(1:121, [idx idx_max]), ll(setdiff(1:121, [idx idx_max])), '^',...
    'MarkerSize', 5, 'MarkerFaceColor', [0    0.4470    0.7410]);
hold on
fig_true = plot(idx, val, 'p','MarkerFaceColor',[0.8500    0.3250    0.0980], 'MarkerSize',16)
hold on
fig_max = plot(idx_max, val_max, 's','MarkerFaceColor', [0.9290    0.6940    0.1250], 'MarkerSize',10)
st = num2str(vec([id0, idx_max]));
dy = 0.05; 
% text(vec(ll_ind(1:20))-0.01, ll(ll_ind(1:20))+dy, st);
text(vec([id0, idx_max])-0.01, ll(vec([id0, idx_max]))+dy, st);

xlim([0 122])
xlabel('indices of candidate Q-matrices')
ylabel('log-likelihood values');
set(gca,'FontSize',12)
pbaspect([16 9 1]);


if idx ~= idx_max
    legend([fig_true, fig_max], 'True Q', 'Alternative Q')
else
    legend([fig_true, fig_max], 'True Q', 'True Q')
end

% print('-r400', 'exhaus_gid_Q5', '-dpng');
% 
% save('gid_Q5.mat')


fprintf('print the true Q generating the data:\n')
Q_arr(:,:,id0)
fprintf('print the estimated Q with the largest log-likelihood:\n')
Q_arr(:,:,idx_max)

ll(id0)-ll(idx_max)


filename = strcat('exhaus_gdina_Q', num2str(id0), '.mat');
save(filename)


% [c_true, c_arr(:, :, idx_max)]
% [g_true, g_arr(:, :, idx_max)]
% [p_true', nu_arr(:, :, idx_max)]


%idvec = get_maxid(Q_arr, 1000, 5, 5)


% K = 2;
% N = 1000;
% 
% Q = [1,0; 0,1; 1,1; 1,1; 1,0];
% J = size(Q, 1);
% 
% % p = [0.1, 0.2, 0.3, 0.4];
% % c = 0.8 * ones(J, 1);
% % g = 0.2 * ones(J, 1);
% 
% pn = rand(1, 4);
% p = pn / sum(pn);
% c = 0.8 * rand(J, 1);
% g = 0.2 * rand(J, 1);
% 
% [X, A] = generate_X(N, p, Q, c, g);
% 
% [nu0, c0, g0, loglik0] = get_cg(X, Q);
% 
% Q_aa = load('Q_aa.mat');
% Q_arr = cell2mat(struct2cell(Q_aa));
% 
% nu_arr = nu0;
% c_arr = c0;
% g_arr = g0;
% ll_arr = zeros(size(Q_arr, 3), 1);
% for r = 1:size(Q_arr, 3)
%     [nu, c, g, loglik] = get_cg(X, Q_arr(:, :, r));
%     nu_arr(:, :, r) = nu;
%     c_arr(:, :, r) = c;
%     g_arr(:, :, r) = g;
%     ll_arr(r) = loglik;
%     fprintf('Just finished iteration %d,\t with loglik %1.4f,\n', ...
%         r, loglik);
% end
% 
% [val, idx] = max(ll_arr);
% 
% % hist(ll_arr)
% % scatter(1:121,ll_arr, 'filled')
% 
% vec = (1:121)';
% 
% 
% %figure(1)
% %fig = figure(1);
% plot(1:121, ll_arr, '^', 'MarkerSize', 5, 'MarkerFaceColor', 'blue');
% %set(gcf, 'Position', [0 0 500 200])
% hold on
% plot(idx, val, 'p','MarkerFaceColor','red', 'MarkerSize',16)
% [ll_val, ll_ind] = sort(ll_arr, 'descend');
% top20 = ll_ind(1:20);
% b = num2str(vec(ll_ind(1:20))); c = cellstr(b);
% dy = 0.5; % displacement so the text does not overlay the data points
% text(vec(ll_ind(1:20))-0.01, ll_arr(ll_ind(1:20))+dy, c);
% hold off
% 
% for i = 1:20
%     Q_arr(:,:,top20(i))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


