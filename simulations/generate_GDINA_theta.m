function [theta, delta] = generate_GDINA_theta(Q)
% This function randomly generate GDINA model parameters given a Q-matrix

[J, K] = size(Q);

S = sum(Q, 2);

p1 = 0.15 + 0.1 * rand(1); 
p2 = 0.85 - 0.1 * rand(1);
delta = [p1*ones(J,1), zeros(J,2^K-1)];
theta = zeros(J, 2^K);


attr_combo = get_I(K, K); % (2^K-1) * K
expand_Q = prod(bsxfun(@power, reshape(Q, [J 1 K]), reshape(attr_combo, [1 2^K-1 K])), 3);
expand_Q = [ones(J,1), expand_Q];
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
    delta(j, ok) = delta_temp(1:2^S(j));

    Aj = [zeros(1, S(j)); get_I(S(j), S(j))];
    % expand Aj to get Mj
    Mj = prod(bsxfun(@power, reshape(Aj, [2^S(j) 1 S(j)]), reshape(Aj, [1 2^S(j) S(j)])), 3);
    theta(j, ok) = (Mj * delta(j, ok)')';
end


end
