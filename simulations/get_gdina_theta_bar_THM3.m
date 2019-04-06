function [Theta_bar, theta_bar_combo, delta_bar_combo, nu_bar, is_monotone] = ...
    get_gdina_theta_bar_THM3(Theta_true, nu_true, Q_bar, expand_Qbar, ind_rng)

% This function constructs alternative model parameters based on Theorem 3
% that are not distinguishable from the true model


rng(ind_rng)

[~, num_class] = size(Theta_true); 
K = log(num_class)/log(2);

attr_combo = get_I(K, K); % (2^K-1) * K

% expand_Qbar = [ones(J,1), prod(bsxfun(@power, reshape(Q_bar, [J 1 K]), ...
%     reshape(attr_combo, [1 2^K-1 K])), 3)];

Theta_bar = Theta_true; nu_bar = zeros(1, num_class);


% rng(1); rng(30); rng(100); rng(1300); rng(1400); rng(1700)
% The following line is the original version
% Theta_bar(1:2, 1:4) = Theta_true(1:2, 1:4) + (-0.1 + 0.2 * rand(2, 4));
% The following is the new version for monotone constraint under Q_bar
% nclass = 2^K; half_class = 2^(K-1);

half1 = 1:(2^(K-1));
half2 = (2^(K-1)+1):(2^K);

A = binary(0:(2^K-1), K);
perturb = zeros(2, 2^(K-1));
A_upright = A(half1, 2:end);
% consider the case with K=5
for kk = 1:(K-1)
    ind_kk = find(sum(A_upright,2)==kk);
    perturb(1:2, ind_kk') = -0.02 + 0.01*(kk-1) + 0.01*rand(2,length(ind_kk));
end
% [-0.15+0.1*rand(2,1), -0.1+0.1*rand(2,1), -0.1+0.1*rand(2,1), -0.05+0.1*rand(2,1)]
Theta_bar(1:2, half1) = Theta_true(1:2, half1) + perturb;


Theta_bar(1, half2) = Theta_true(1, half1) + ...
    (Theta_true(1, half2) - Theta_true(1, half1)) .* ...
    (Theta_true(2, half2) - Theta_bar(2,half1)) .* (nu_true(half2))' ./ ...
    ((Theta_true(2,half1) - Theta_bar(2,half1)) .* (nu_true(half1))' + ...
    (Theta_true(2, half2) - Theta_bar(2,half1)) .* (nu_true(half2))');

Theta_bar(2, half2) = Theta_true(2, half1) + ...
    (Theta_true(2, half2) - Theta_true(2, half1)) .* ...
    (Theta_true(1,half2) - Theta_bar(1,half1)) .* (nu_true(half2))' ./ ...
    ((Theta_true(1,half1) - Theta_bar(1,half1)) .* (nu_true(half1))' + ...
    (Theta_true(1,half2) - Theta_bar(1,half1)) .* (nu_true(half2))');

nu_bar(half2) = ( (Theta_true(2,half1) - Theta_bar(2,half1)) .* ...
    nu_true(half1)' + (Theta_true(2,half2) - Theta_bar(2,half1)) .* nu_true(half2)' )  ./ ...
    (Theta_bar(2,half2) - Theta_bar(2,half1));


nu_bar(half1) = nu_true(half1)' + nu_true(half2)' - nu_bar(half2);
nu_bar = nu_bar';


% alternative parameters under Q_bar
theta_bar_combo = Theta_bar(:, [1; bin2ind(attr_combo)+1]) .* expand_Qbar;
% check if alternative parameters respect monotone constraint under Q_bar
[is_monotone, delta_bar_combo] = check_monotone(Q_bar, theta_bar_combo);



end