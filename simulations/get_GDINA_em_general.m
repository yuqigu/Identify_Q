function [theta_sub, delta_sub, nu, loglik] = get_GDINA_em_general(X, Q, A, ...
    expand_Q_sub, attr_patt_eff, theta_ini_sub, expand_S)

% EM algorithm for estimating category proportions nu_c and
% GDINA parameters theta, delta

% @param X               : N x J response matrix
% @param Q               : J x K Q-matrix
% @param theta_ini       : initilization of GDINA paramters theta
% @param nu_ini          : initilization of proportions




[J, K] = size(Q);

N = size(X, 1);
thres = 0.5/N;

resp_vecs = X;
count_resp = ones(N, 1);
prop_resp = count_resp/N;


C_eff = size(attr_patt_eff,1);
theta_sub = theta_ini_sub;


% 2^K * 2^K, rows in bin-vector order, columns in one-order
% expand_A = [ones(size(A,1), 1), prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), ...
%     reshape(attr_combo, [1 2^K-1 K])), 3)];

% size 2^K * size(attr_patt_eff,1)
expand_A_sub = prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), ...
    reshape(attr_patt_eff, [1 C_eff K])), 3);
C = size(A, 1);
% pa_nume = ones(1, C)/C;
% always initialize uniformly
pa_nume = ones(1, C);

ok_ex_all = (expand_Q_sub==1); % logical version of expand_Q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_sub = zeros(J, C_eff);

S = sum(Q, 2); % row sums of Q

MM = cell(J, 1);
inv_MM = cell(J, 1);
for j = 1:J
    MM{j} = prod(bsxfun(@power, reshape(expand_S{j}, ...
        [2^S(j) 1 S(j)]), reshape(expand_S{j}, [1 2^S(j) S(j)])), 3);
    inv_MM{j} = inv(MM{j});
end

for j = 1:J
    delta_sub(j, ok_ex_all(j, :)) = ( inv_MM{j} * theta_sub(j, ok_ex_all(j, :))' )';
end

% positive response probability of each attribute pattern to each item, J*(2^K)
% columns arranged in the binary vector order
prob_posi = delta_sub * expand_A_sub'; 

itera = 0;

err = 1;


while err > 1e-4 && itera < 1000
    old_theta_sub = theta_sub;

    posi_part = prod( bsxfun(@power, prob_posi', reshape(resp_vecs', [1 J N])), 2 );
    nega_part = prod( bsxfun(@power, 1 - prob_posi', reshape(1 - resp_vecs', [1 J N])), 2 );

    alpha_li0 = (squeeze(posi_part)' .* squeeze(nega_part)') .* pa_nume; % unnormalized version phi
    
    alpha_li = alpha_li0 ./ sum(alpha_li0, 2); % size N * 2^K    
    
    pa_nume = sum(N * prop_resp .* alpha_li, 1);
    % pa_nume = lambda + numer;
    % dprior = real(lambda + pow * numer);
    
    for j = 1:J
        Aj = expand_S{j};
        Rj = zeros(2^S(j), 1);
        Ij = zeros(2^S(j), 1);
        for jj = 1:2^S(j)
            ind_jj = prod( A(:, Q(j,:)==1) == repmat(Aj(jj,:), [C,1]), 2 );
            Rj(jj) = sum(alpha_li(:, ind_jj==1)' * (count_resp .* resp_vecs(:, j)));
            Ij(jj) = sum(alpha_li(:, ind_jj==1)' * count_resp);
        end
        
        ratio = zeros(2^S(j), 1);
        ratio(Ij > 0) = Rj(Ij > 0) ./ Ij(Ij > 0);
        theta_sub(j, ok_ex_all(j, :)) = ratio;

        
        delta_sub(j, ok_ex_all(j, :)) = ( inv_MM{j} * theta_sub(j, ok_ex_all(j, :))' )';
    end
    prob_posi = delta_sub * expand_A_sub';
    
    itera = itera + 1;
    nu = pa_nume' / sum(pa_nume);
    nu = nu / sum(nu);
    
    err = max(max(abs(old_theta_sub - theta_sub)));
    fprintf('EM: iteration %d,\t size above threshold %d,\t with Err %1.6f \n',...
        itera, sum(nu>thres), err);
end


loglik = N*sum(prop_resp .* log(prod(bsxfun(@power, reshape(prob_posi',[1 C J]), ...
    reshape(resp_vecs, [N 1 J])) .* ...
    bsxfun(@power, 1-reshape(prob_posi',[1 C J]), ...
    1-reshape(resp_vecs, [N 1 J])),3)*nu));

end
