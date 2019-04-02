function [theta, delta, p, loglik, is_mono, is_mono_str] = get_GDINAprob_mono(X, Q, theta_ini, p_ini, expand_S)

% This function uses the EM algorithm to compute the MLE under the GDINA
% model; and in the last step further checks whether the estimated parameters 
% satisfy the monotonicity constraint.

% @param X               : N x J response matrix
% @param Q               : J x K Q-matrix
% @param theta_ini       : initilization of GDINA paramters theta
% @param p_ini          : initilization of proportions

% @return nu_c : catogory proportions
% @return delta: GDINA paramters

% EM algorithm for estimating category proportions nu_c and
% GDINA parameters theta, delta

[J, K] = size(Q);
N = size(X, 1);

resp_vecs = X;
count_resp = ones(N, 1);
prop_resp = count_resp/N;

theta = theta_ini;
p = p_ini';
old_theta = ones(J, 2^K);


attr_combo = get_I(K, K); % (2^K-1) * K
expand_Q = [ones(J,1), prod(bsxfun(@power, reshape(Q, [J 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];
A = binary(0:(2^K-1), K);
expand_A = [ones(size(A,1), 1), prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];

ok_ex_all = (expand_Q==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = zeros(J, 2^K);

S = sum(Q, 2);
MM = cell(J, 1);
inv_MM = cell(J, 1);
for j = 1:J
    MM{j} = prod(bsxfun(@power, reshape(expand_S{j}, ...
        [2^S(j) 1 S(j)]), reshape(expand_S{j}, [1 2^S(j) S(j)])), 3);
    inv_MM{j} = inv(MM{j});
end

for j = 1:J
    delta(j, ok_ex_all(j, :)) = ( inv_MM{j} * theta(j, ok_ex_all(j, :))' )';
end
prob_posi = delta * expand_A';

max_itera = 5000;
itera = 0;

while sum(sum(abs(old_theta - theta))) > 1e-5 && itera<max_itera
    old_theta = theta;

    posi_part = prod( bsxfun(@power, prob_posi', reshape(resp_vecs', [1 J N])), 2 );
    nega_part = prod( bsxfun(@power, 1 - prob_posi', reshape(1 - resp_vecs', [1 J N])), 2 );
    
    alpha_li0 = squeeze(posi_part)' .* squeeze(nega_part)' .* repmat(p', [N, 1]);
    
    alpha_li = alpha_li0 ./ repmat(sum(alpha_li0, 2), [1, 2^K]);
    
    
    numer = sum(repmat(count_resp, [1, 2^K]) .* alpha_li, 1);
    p = numer / sum(numer);
    p = p';
    
    for j = 1:J
        Aj = expand_S{j};
        Rj = zeros(2^S(j), 1);
        Ij = zeros(2^S(j), 1);
        for jj = 1:2^S(j)
            ind_jj = prod( A(:, Q(j,:)==1) == repmat(Aj(jj,:), [2^K,1]), 2 );
            
            Rj(jj) = sum(alpha_li(:, ind_jj==1)' * (count_resp .* resp_vecs(:, j)));
            Ij(jj) = sum(alpha_li(:, ind_jj==1)' * count_resp);
        end
        
        ratio = zeros(2^S(j), 1);
        ratio(Ij > 0) = Rj(Ij > 0) ./ Ij(Ij > 0);
        theta(j, ok_ex_all(j, :)) = ratio;

        
        delta(j, ok_ex_all(j, :)) = ( inv_MM{j} * theta(j, ok_ex_all(j, :))' )';
    end
    prob_posi = delta * expand_A';
    
    itera = itera + 1;
    % fprintf('Finished EM iteration %d, with error %1.6f\n', itera, sum(sum(abs(old_theta - theta))));
%     if itera == 370
%         break
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following code checks if theta satisfies the monotonicity constraint
% specified by Q-matrix

Gamma = get_ideal_resp(Q, A);

% obtain the posi. resp. probability of the all-one attribute pattern to
% all the items; this should
% The posi. resp. prob of capable patterns for each item
capable_prob = prob_posi(:,end);
% The highest posi. resp. prob of non-capable patterns for each item
nonca_max_prob = zeros(J, 1); 
for j = 1:J
    % capable_prob(j) = max(prob_posi(j, Gamma(j,:)==1));
    nonca_max_prob(j) = max(prob_posi(j, Gamma(j,:)==0));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% officially check if monotonicity constraint is satisfied
is_mono = all(capable_prob > nonca_max_prob);

% a further decision about monotonicity (may be stronger than required)
% check if every nonzero element in delta is positive (can't be negative)
% is_mono_str = all(all(delta >= 0));

is_mono_str = all(delta(delta ~= 0) > 0);


loglik = N*sum(prop_resp.*log(prod(bsxfun(@power, reshape(prob_posi',[1 2^K J]), ...
    reshape(resp_vecs, [N 1 J])) .* ...
    bsxfun(@power, 1-reshape(prob_posi',[1 2^K J]), ...
    1-reshape(resp_vecs, [N 1 J])),3)*p));

end
