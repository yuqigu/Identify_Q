function [nu_c, c_i, g_i, loglik] = get_cg(X, Q, p_ini, c_ini, g_ini)

% @param X               : n x m response matrix
% @param Q               : m x k Q-matrix
%
% @return c_i            : m x 1 slipping vector
% @return g_i            : m x 1 slipping vector

% EM algorithm for estimating c, g and mixing coeffcients \pi

[M, K] = size(Q);


N = size(X, 1);
resp_vecs = unique(X, 'rows');
prop_resp = histc(bin2ind(X)+1, 1:2^M)/N; % 1:2^M is binranges
prop_resp = prop_resp(prop_resp>0); % proportion of indices of subsets of {all possible items combinations} in X
N_x = size(resp_vecs,1); % number of unique subsets of items in X
    

all_attr = reshape(binary(0:(2^K-1),K),[1 2^K K]);
Q = reshape(Q, [M 1 K]);
xi_ic = prod(bsxfun(@power, all_attr, Q), 3); % ideal response matrix

old_c = ones(M,1);
old_g = zeros(M,1);
nu_c = ones(1,1,2^K)/2^K; % initialization of mixing coefficients (profile proportions)
nu_c(:,:,1:2^K) = p_ini;
c_i = c_ini; % initialization of 1-slipping parameter
g_i = g_ini; % initialization of guessing parameter


while (max(abs(old_c-c_i))+max(abs(old_g-g_i)))>1e-6
    
    pi_ic = bsxfun(@times, c_i, xi_ic) + bsxfun(@times, g_i, 1-xi_ic); % M * 2^K, prob matrix for positive responses
    
    % prob for each response pattern in X and each attribute profile, under
    % current values of g and c
    alpha_rc = bsxfun(   @times, reshape(nu_c,[1 1 2^K]), ...
        prod(  bsxfun(@power, reshape(pi_ic, [1 M 2^K]), resp_vecs) .* ...
        bsxfun(@power, 1-reshape(pi_ic,[1 M 2^K]), 1-resp_vecs), 2  )   ); 
    
    % N_x * 1 * 2^K, prob for each response pattern in X and each attribute profile
    alpha_rc = bsxfun(@rdivide, alpha_rc, sum(alpha_rc, 3)); 
    % alpha_rc = alpha_rc ./ sum(alpha_rc, 3);
    
    nu_c = sum(bsxfun(@times, alpha_rc, prop_resp), 1);
    
    % N_x * M, marginal prob for each response pattern 
    prob_know_ri = sum(bsxfun(@times, alpha_rc, ...
        reshape(xi_ic, [1 M 2^K])), 3);
    
    % prob_know_ri0 = reshape(alpha_rc, [N_x 2^K]) * reshape(xi_ic, [2^K M]);
    
    old_c = c_i;
    old_g = g_i;
    
    c_i = sum(bsxfun(@times,resp_vecs.*prob_know_ri,prop_resp),1)'./...
        sum(bsxfun(@times,prob_know_ri,prop_resp),1)';
    
    g_i = sum(bsxfun(@times,resp_vecs.*(1-prob_know_ri),prop_resp),1)'./...
        sum(bsxfun(@times,1-prob_know_ri,prop_resp),1)';
end

pi_ic = bsxfun(@times, c_i, xi_ic) + bsxfun(@times, g_i, 1-xi_ic);

alpha_rc = bsxfun(@times, reshape(nu_c,[1 1 2^K]), ...
    prod(bsxfun(@power, reshape(pi_ic,[1 M 2^K]), resp_vecs) .* ...
    bsxfun(@power, 1-reshape(pi_ic,[1 M 2^K]), 1-resp_vecs),2));

alpha_rc = bsxfun(@rdivide, alpha_rc, sum(alpha_rc,3));

nu_c = sum(bsxfun(@times,alpha_rc,prop_resp), 1);
nu_c = reshape(nu_c, [2^K 1]);

loglik = N*sum(prop_resp.*log(prod(bsxfun(@power, reshape(pi_ic',[1 2^K M]), ...
    reshape(resp_vecs, [N_x 1 M])) .* ...
    bsxfun(@power, 1-reshape(pi_ic',[1 2^K M]), ...
    1-reshape(resp_vecs, [N_x 1 M])),3)*nu_c));

end