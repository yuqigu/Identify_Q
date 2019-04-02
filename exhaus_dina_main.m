%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DINA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_aa = load('Q_aa.mat');
Q_arr = cell2mat(struct2cell(Q_aa));

% vec = (1:121)';
% ind_id = [15, 17, 18, 23, 24, 26, 33, 35, 36, 39,...
%     43, 45, 47, 48, 49, 51, 52, 53, 59, 60, 62, 65, 66,...
%     67, 69, 70, 71, 74, 76, 77, 86, 87, 89, 92, 93, 94, ...
%     96, 97, 98, 101, 103, 104, 110, 112, 113]';
% 
% 
% ind_not_lid = [95];


%%% the following are several simulation settings for the paper %%%
% simulation setting 1: ABC satisfied; two cases
% strictly identifiable
ind_strid_two_I2 = [15, 17, 23];
ind_strid_one_I2 = [18, 24, 26];

% simulation setting 2: AB satisfied, Thm 2(b); only one case
% generically identifiable (globally)
ind_gid = [5, 11, 13];

% simulation setting 3: attr>=2 and complete, but v=1 in Theorem 2(a)
% not even locally generically identifiable
ind_nlid_thm2a = [6, 8, 16];

% simulation setting 4: attr>=2 but not complete
% not even locally generically identifiable
id_nlid_incom = [9, 21, 25];

% simulation setting 5: complete but attr=1
% not even locally generically identifiable
id_nlid_attr1 = [2, 4, 10];


% simulation setting 6: complete but attr=1
% not even locally generically identifiable
id_nlid_incom_attr1 = [3, 19, 55];




clear;
load('exhaus_dina_Q55.mat')

val = ll(id0);
[val_max, idx_max] = max(ll);

figure
fig_all = plot(setdiff(1:121, [id0 idx_max]), ll(setdiff(1:121, [id0 idx_max])), '^',...
    'MarkerSize', 5, 'MarkerFaceColor', [0    0.4470    0.7410]);
hold on
fig_max = plot(idx_max, val_max, 's','MarkerFaceColor', [0.4940    0.1840    0.5560],...
    'MarkerEdgeColor', [0.4940    0.1840    0.5560], 'MarkerSize',15);
hold on
fig_true = plot(id0, val, 'p','MarkerFaceColor',[0.8500    0.3250    0.0980], ...
    'MarkerEdgeColor', [0.8500    0.3250    0.0980], 'MarkerSize',20);
st = num2str(vec([id0, idx_max]));
text(vec([id0, idx_max]), ll(vec([id0, idx_max]))+0.05, st);
xlim([0 122])
xlabel('indices of candidate Q-matrices')
ylabel('log-likelihood values');
set(gca,'FontSize',12)
pbaspect([16 9 1]);
legend([fig_true, fig_max], 'true Q', 'the Q with largest likelihood',...
    'location', 'northeast')


print('-r400', 'exhaus_dina_Q55', '-dpng');




% if idx ~= idx_max
%     legend([fig_true, fig_max], 'True Q', 'Alternative Q')
% else
%     legend([fig_true, fig_max], 'True Q', 'True Q')
% end


save('gid_Q5.mat')


fprintf('print the true Q generating the data:\n')
Q_arr(:,:,id0)
fprintf('print the estimated Q with the largest log-likelihood:\n')
Q_arr(:,:,idx_max)

ll(id0)-ll(idx_max)


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

% % the index of the true Q-matrix used to generate the data
% id0 = 5;
% 
% 
% 
% 
% 
% Q0 = Q_arr(:,:,id0);
% 
% [J, K] = size(Q0);
% 
% % p = [0.1, 0.2, 0.3, 0.4];
% % c = 0.8 * ones(J, 1);
% % g = 0.2 * ones(J, 1);
% 
% rng(0)
% pn = gamrnd(5 * ones(1,2^K), ones(1,2^K));
% % [rand(1), rand(1)*4, rand(1)*3, rand(1)*2];
% p_true = pn / sum(pn);
% c_true = 0.9 - 0.2 * rand(J, 1);
% g_true = 0.1 + 0.2 * rand(J, 1);
% 
% % set initial values
% p_ini = p_true; c_ini = c_true; g_ini = g_true;
% 
% 
% N = 10000;
% [X, A] = generate_X(N, p_true, Q0, c_true, g_true);
% % [nu0, c0, g0, loglik0] = get_cg(X, Q0, p_ini, c_ini, g_ini);
% 
% % set the number of random initializations
% num_ran = 5;
% 
% 
% nu_arr = p_true;
% c_arr = c_true;
% g_arr = g_true;
% Q = Q0;
% 
% ll = zeros(size(Q_arr, 3), 1);
% for r = 1:size(Q_arr, 3)
%     
%     % Revision: 5 random initializations
%     llk_vec = zeros(num_ran, 1);
%     p_vec = zeros(2^K, num_ran);
%     c_vec = zeros(J, num_ran);
%     g_vec = zeros(J, num_ran);
% 
% %     [p_hat0, c_hat0, g_hat0, loglike0] = get_cg(X, Q_arr(:,:,r),...
% %             p_true, c_true, c_true);
% %         
% %     p_vec(:,1) = p_hat0; c_vec(:,1)=c_hat0; g_vec(:,1)=g_hat0;
% %     llk_vec(1) = loglike0;
%         
%     for ii = 1:num_ran
%         pn0 = gamrnd(5 * ones(1,2^K), ones(1,2^K));
%         p_ini = pn0 / sum(pn0);
%         c_ini = 0.85 - 0.1 * rand(J, 1);
%         g_ini = 0.15 + 0.1 * rand(J, 1);
% 
%         [p_hat, c_hat, g_hat, loglike] = get_cg(X, Q_arr(:,:,r),...
%             p_ini, c_ini, g_ini);
% 
%         p_vec(:,ii) = p_hat; c_vec(:,ii)=c_hat; g_vec(:,ii)=g_hat;
% 
%         llk_vec(ii) = loglike;
%     end
%     [~, ind] = max(llk_vec);
%     
%     % [nu, c, g, loglik] = get_cg(X, Q_arr(:, :, r), p_ini, c_ini, g_ini);
%     nu_arr(:, :, r) = p_vec(:,ind);
%     c_arr(:, :, r) = c_vec(:,ind);
%     g_arr(:, :, r) = g_vec(:,ind);
%     ll(r) = llk_vec(ind);
%     
%     fprintf('Just finished iteration %d,\t with loglik %1.4f,\n', ...
%         r, llk_vec(ind));
% end
