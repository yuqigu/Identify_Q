function [theta, delta, p, loglik] = get_GDINAprob(X, Q, theta_ini, p_ini, expand_S)
% function [nu, delta] = get_GDINAprob(X, Q, theta_ini)

% @param X               : N x J response matrix
% @param Q               : J x K Q-matrix
% @param theta_ini       : initilization of GDINA paramters theta
% @param p_ini          : initilization of proportions

% @return nu_c : catogory proportions
% @return delta: GDINA paramters

% EM algorithm for estimating category proportions nu_c and
% GDINA parameters theta, delta

[M, K] = size(Q);
N = size(X, 1);

resp_vecs = X;
count_resp = ones(N, 1);
prop_resp = count_resp/N;

theta = theta_ini;
p = p_ini';
old_theta = ones(M, 2^K);


attr_combo = get_I(K, K); % (2^K-1) * K
expand_Q = [ones(M,1), prod(bsxfun(@power, reshape(Q, [M 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];
A = binary(0:(2^K-1), K);
expand_A = [ones(size(A,1), 1), prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];

ok_ex_all = (expand_Q==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = zeros(M, 2^K);

S = sum(Q, 2);
MM = cell(M, 1);
inv_MM = cell(M, 1);
for j = 1:M
    MM{j} = prod(bsxfun(@power, reshape(expand_S{j}, ...
        [2^S(j) 1 S(j)]), reshape(expand_S{j}, [1 2^S(j) S(j)])), 3);
    inv_MM{j} = inv(MM{j});
end

for j = 1:M
    delta(j, ok_ex_all(j, :)) = ( inv_MM{j} * theta(j, ok_ex_all(j, :))' )';
end
prob_posi = delta * expand_A';

max_itera = 5000;
itera = 0;

while sum(sum(abs(old_theta - theta))) > 1e-5 && itera<max_itera
    old_theta = theta;

    posi_part = prod( bsxfun(@power, prob_posi', reshape(resp_vecs', [1 M N])), 2 );
    nega_part = prod( bsxfun(@power, 1 - prob_posi', reshape(1 - resp_vecs', [1 M N])), 2 );
    
    % alpha_li0 = squeeze(posi_part)' .* squeeze(nega_part)' .* p';
    alpha_li0 = squeeze(posi_part)' .* squeeze(nega_part)' .* repmat(p', [N, 1]);
    
    % alpha_li = alpha_li0 ./ sum(alpha_li0, 2);
    alpha_li = alpha_li0 ./ repmat(sum(alpha_li0, 2), [1, 2^K]);
    
    % p = sum(count_resp .* alpha_li, 1) ./ sum(sum(count_resp .* alpha_li, 1)); % 2^K * 1
    % p = p';
    
    numer = sum(repmat(count_resp, [1, 2^K]) .* alpha_li, 1);
    p = numer / sum(numer);
    p = p';
    
    for j = 1:M
        % Q_j_try = [zeros(1, S(j)); get_I(S(j), S(j))];
        Aj = expand_S{j};
        % ok_ex = (expand_Q(j, :)==1);
        % ok_ex = ok_ex_all(j, :);
        Rj = zeros(2^S(j), 1);
        Ij = zeros(2^S(j), 1);
        for jj = 1:2^S(j)
            % 0603 Saturday change
            % ind_jj = prod( A(:, Q(j,:)==1) == Aj(jj,:), 2 );
            ind_jj = prod( A(:, Q(j,:)==1) == repmat(Aj(jj,:), [2^K,1]), 2 );
            
            Rj(jj) = sum(alpha_li(:, ind_jj==1)' * (count_resp .* resp_vecs(:, j)));
            Ij(jj) = sum(alpha_li(:, ind_jj==1)' * count_resp);
        end
        
        %%%% 0613 change starts %%%%%
        ratio = zeros(2^S(j), 1);
        ratio(Ij > 0) = Rj(Ij > 0) ./ Ij(Ij > 0);
        theta(j, ok_ex_all(j, :)) = ratio;
        %%%% 0613 change ends %%%%%
                
        % Aj = [zeros(1, S(j)); get_I(S(j), S(j))];
        % expand Aj to get Mj
        % Mj = prod(bsxfun(@power, reshape(Aj, [2^S(j) 1 S(j)]), reshape(Aj, [1 2^S(j) S(j)])), 3);
        % delta(j, ok_ex_all(j, :)) = ( Mj \ theta(j, ok_ex_all(j, :))' )';
        
        delta(j, ok_ex_all(j, :)) = ( inv_MM{j} * theta(j, ok_ex_all(j, :))' )';
    end
    prob_posi = delta * expand_A';
    
    itera = itera + 1;
    % fprintf('Finished EM iteration %d, with error %1.6f\n', itera, sum(sum(abs(old_theta - theta))));
%     if itera == 370
%         break
%     end
end

% p = reshape(p, [2^K 1]);

loglik = N*sum(prop_resp.*log(prod(bsxfun(@power, reshape(prob_posi',[1 2^K M]), ...
    reshape(resp_vecs, [N 1 M])) .* ...
    bsxfun(@power, 1-reshape(prob_posi',[1 2^K M]), ...
    1-reshape(resp_vecs, [N 1 M])),3)*p));

end
