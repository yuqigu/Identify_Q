function [Rp] = get_rp_GDINA(Q, I_full, prob_posi)
% produce R-matrix for all rows denoted by I, storing proportions of 
% positive response for each response pattern in I

[M, K] = size(Q);
% N_I = 2^M-1; % full I

I_full1 = [zeros(1,M); I_full]; % (n + 1) * M

% B1 = zeros(M,size(A,1));
% for i=1:M
%     for j  = 1:size(A,1);
%         B1(i,j) = prod(A(j,:).^Q(i,:));
%     end
% end

% %  M * 2^K, indicators of capabilities, ideal responses for all items and all attribute profiles
% B1 = prod( bsxfun(@power, reshape(A, [1 size(A,1) K]), reshape(Q, [M 1 K])), 3);
% 
% Bs_r = bsxfun(@plus, g, bsxfun(@times, B1, c-g)); % M * 2^K prob of positive response for the capable
% Bs_w = 1-Bs_r; % M * 2^K prob of negative response for the capable

% attr_combo = get_I(K, K); % (2^K-1) * K
% expand_A = prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), reshape(attr_combo, [1 2^K-1 K])), 3); 
% expand_A = [ones(size(A,1), 1), expand_A]; % M matrix in the GDINA paper
% % prob of positive response for all single items and all attribute profiles, M * 2^K
% prob_posi = delta * expand_A';
prob_nega = 1 - prob_posi;

% produce R-matrix for all rows denoted by I, storing proportions of
% positive response for each response pattern in I
%%%%% arrayfun is slower than for loop %%%%%
% Rp0 = cell2mat(arrayfun(@(i) prod(Bs_w(I_full1(i,:)==0,:),1) .* prod(Bs_r(I_full1(i,:)==1,:),1),...
%         (1:size(I_full1,1))', 'UniformOutput', false));

Rp = zeros(2^M, 2^K);
for i = 1:size(I_full1, 1)
    Rp(i, :) = prod(prob_nega(I_full1(i,:)==0,:),1) .* prod(prob_posi(I_full1(i,:)==1,:),1);
end


end