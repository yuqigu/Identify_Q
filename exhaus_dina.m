function [] = exhaus_dina(id0, N)


% the index of the true Q-matrix used to generate the data
% id0 = 99; N=10^5;

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
% 
% ind_gid = [5, 11];


Q0 = Q_arr(:,:,id0);

[J, K] = size(Q0);


rng(0)
pn = gamrnd(5 * ones(1,2^K), ones(1,2^K));
p_true = pn / sum(pn);
c_true = 0.9 - 0.2 * rand(J, 1);
g_true = 0.1 + 0.2 * rand(J, 1);

% set initial values
% p_ini = p_true; c_ini = c_true; g_ini = g_true;


[X, ~] = generate_X(N, p_true, Q0, c_true, g_true);

% set the number of random initializations
num_ran = 5;


nu_arr = p_true;
c_arr = c_true;
g_arr = g_true;
% Q = Q0;

ll = zeros(size(Q_arr, 3), 1);
for r = 1:size(Q_arr, 3)
    
    % Revision: 5 random initializations
    llk_vec = zeros(num_ran, 1);
    p_vec = zeros(2^K, num_ran);
    c_vec = zeros(J, num_ran);
    g_vec = zeros(J, num_ran);

    for ii = 1:num_ran
        pn0 = gamrnd(5 * ones(1,2^K), ones(1,2^K));
        p_ini = pn0 / sum(pn0);
        c_ini = 0.85 - 0.1 * rand(J, 1);
        g_ini = 0.15 + 0.1 * rand(J, 1);

        [p_hat, c_hat, g_hat, loglike] = get_cg(X, Q_arr(:,:,r),...
            p_ini, c_ini, g_ini);

        p_vec(:,ii) = p_hat; c_vec(:,ii)=c_hat; g_vec(:,ii)=g_hat;

        llk_vec(ii) = loglike;
    end
    [~, ind] = max(llk_vec);
    
    nu_arr(:, :, r) = p_vec(:,ind);
    c_arr(:, :, r) = c_vec(:,ind);
    g_arr(:, :, r) = g_vec(:,ind);
    ll(r) = llk_vec(ind);
    
    fprintf('Just finished iteration %d,\t with loglik %1.4f,\n', ...
        r, llk_vec(ind));
end

filename = strcat('exhaus_dina_Q', num2str(id0), '_N', num2str(N), '.mat');
save(filename)


end