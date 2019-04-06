% This file performs exhaustive search of Q-matrix in the set of 5*2
% candidate Q-matrices under GDINA model, corresponding to simulation studies 
% V and VI in the Supplementary Material

load('exhaus_gdina_Q9.mat')

val = ll(id0);
vec = (1:size(Q_arr,3))';

ind_strmono = vec(is_monostr_arr==1);
ind_NONmono = setdiff(vec, ind_strmono);

ll_non0_onlymono = ll;
ll_non0_onlymono(ind_NONmono) = -Inf;
[val_max, idx_max] = max(ll_non0_onlymono);

figure
fig_all = plot(ind_strmono, ll(ind_strmono), '^', ...
    'MarkerSize', 10, 'MarkerFaceColor', [0    0.4470    0.7410]);
hold on
fig_max = plot(idx_max, val_max, 's','MarkerFaceColor', ...
    [0.4940    0.1840    0.5560], 'MarkerEdgeColor', ...
    [0.4940    0.1840    0.5560], 'MarkerSize', 15);
hold on
fig_true = plot(id0, val, 'p','MarkerFaceColor', ...
    [0.8500    0.3250    0.0980], 'MarkerEdgeColor', ...
    [0.8500    0.3250    0.0980],  'MarkerSize', 20);
st = num2str(vec([id0, idx_max]));
dy = 0.1; 
text(vec([id0, idx_max]), ll(vec([id0, idx_max])), st);
xlim([0 122])
xlabel('indices of candidate Q-matrices')
ylabel('log-likelihood values');
set(gca,'FontSize',12)
pbaspect([16 9 1]);
legend([fig_true, fig_max], 'true Q', 'the Q with largest likelihood',...
    'Location','northeast');
print('-r400', 'exhaus_gdina_onlymono_Q9', '-dpng');



%%%%%%%%%%%% check if Theorem 4 conditions D, E are satisfied %%%%%%%%%%
num_Q = size(Q_arr,3);
cond_DE_vec = zeros(num_Q, 1);
for r = 1:num_Q
    [cond_hold, ~, ~, ~] = check_Theorem4(Q_arr(:,:,r));
    cond_DE_vec(r) = cond_hold;
end





% %%%%%%%%%%% check old versions %%%%%%%%%%%%%%%
% % if idx ~= idx_max
% %     legend([fig_true, fig_max], 'True Q', 'Alternative Q')
% % else
% %     legend([fig_true, fig_max], 'True Q', 'True Q')
% % end
% 
% load('exhaus_gdina_Q15_old.mat')
% Q0
% ll(id0) - ll(idx_max)
% 
% 
% [ll_sort, ind_sort] = sort(ll,'descend');
% Q_arr(:,:,ind_sort(1:15)) - Q_arr(:,:,id0)
% 
% ind_max15_ll = ind_sort(1:15);
% is_mono_arr(ind_max15_ll)
% 
% print('-r400', 'exhaus_gdina_Q15', '-dpng');
% 
% % save('exhaus_gdina_Q5.mat')
% 
% 
% 
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % check monotonicity constraint
% % A = binary(2^(2^K-1), K);
% % is_mono = zeros(size(Q_arr,3), 1);
% % for r = 1:size(Q_arr,3)
% %     Q = Q_arr(:,:,r);
% %     
% %     prob_posi = 
% %     
% %     % get the ideal response matrix Gamma
% %     Gamma = get_ideal_resp(Q, A);
% % 
% %     % obtain the posi. resp. probability of the all-one attribute pattern to
% %     % all the items; this should
% %     % The posi. resp. prob of capable patterns for each item
% %     capable_prob = prob_posi(:,end);
% %     % The highest posi. resp. prob of non-capable patterns for each item
% %     nonca_max_prob = zeros(J, 1); 
% %     for j = 1:J
% %         % capable_prob(j) = max(prob_posi(j, Gamma(j,:)==1));
% %         nonca_max_prob(j) = max(prob_posi(j, Gamma(j,:)==0));
% %     end
% % 
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % officially check if monotonicity constraint is satisfied
% %     is_mono(r) = all(capable_prob > nonca_max_prob);
% % 
% % 
% % end
