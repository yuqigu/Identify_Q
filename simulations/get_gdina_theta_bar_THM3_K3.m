function [Theta_bar, theta_bar_combo, delta_bar_combo, nu_bar, is_monotone] = ...
    get_gdina_theta_bar_THM3_K3(Theta_true, nu_true, Q_bar, expand_Qbar, ind_rng)

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
Theta_bar(1:2, 1:4) = Theta_true(1:2, 1:4) + ...
    [-0.15+0.1*rand(2,1), -0.1+0.1*rand(2,1), -0.1+0.1*rand(2,1), -0.05+0.1*rand(2,1)];


Theta_bar(1, 5:8) = Theta_true(1, 1:4) + (Theta_true(1, 5:8) - Theta_true(1, 1:4)) .* (Theta_true(2,5:8) - Theta_bar(2,1:4)) .* (nu_true(5:8))' ./ ...
    ((Theta_true(2,1:4) - Theta_bar(2,1:4)) .* (nu_true(1:4))' + (Theta_true(2,5:8) - Theta_bar(2,1:4)) .* (nu_true(5:8))');

Theta_bar(2, 5:8) = Theta_true(2, 1:4) + (Theta_true(2, 5:8) - Theta_true(2, 1:4)) .* (Theta_true(1,5:8) - Theta_bar(1,1:4)) .* (nu_true(5:8))' ./ ...
    ((Theta_true(1,1:4) - Theta_bar(1,1:4)) .* (nu_true(1:4))' + (Theta_true(1,5:8) - Theta_bar(1,1:4)) .* (nu_true(5:8))');

nu_bar(5:8) = ( (Theta_true(2,1:4) - Theta_bar(2,1:4)) .* nu_true(1:4)' + (Theta_true(2,5:8) - Theta_bar(2,1:4)) .* nu_true(5:8)' )  ./ ...
    (Theta_bar(2,5:8) - Theta_bar(2,1:4));
nu_bar(1:4) = nu_true(1:4)' + nu_true(5:8)' - nu_bar(5:8);
nu_bar = nu_bar';


% alternative parameters under Q_bar
theta_bar_combo = Theta_bar(:, [1; bin2ind(attr_combo)+1]) .* expand_Qbar;
% check if alternative parameters respect monotone constraint under Q_bar
[is_monotone, delta_bar_combo] = check_monotone(Q_bar, theta_bar_combo);



end