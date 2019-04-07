function [p_mse, c_mse, g_mse] = get_Q24_single_MSE(N, p_true, c_true, g_true)

% This function computes the MLE of DINA model under the 4*2 Q-matrix, for
% a given set of true parameters (p_true, c_true, g_true)
% Revision: use EM with 10 random starts and select the highest MLE

num_ran = 10;


Q = [1 0; 0 1; 1 0; 0 1];

[J, K] = size(Q);


%%% simulations start %%%
Nset = size(p_true, 2); % number of true parameter sets
Niter = 200; % number of iterations per true parameter set


p_mse = zeros(Nset, 1);
c_mse = zeros(Nset, 1);
g_mse = zeros(Nset, 1);

rng('default');
for rr = 1:Nset
    
    p_dist = zeros(Niter, 1);
    g_dist = p_dist; c_dist = p_dist;

    for iter = 1:Niter
        [X, ~] = generate_X(N, p_true(:,rr), Q, c_true(:,rr), g_true(:,rr));
        
        % Revision: 10 random initializations
        llk_vec = zeros(num_ran, 1);
        p_vec = zeros(2^K, num_ran);
        c_vec = zeros(J, num_ran);
        g_vec = zeros(J, num_ran);
        
        for ii = 1:num_ran
                pn0 = gamrnd(3*ones(1,2^K), ones(1,2^K));
                p_ini = pn0 / sum(pn0);
                c_ini = 0.9 - 0.2 * rand(J, 1);
                g_ini = 0.1 + 0.2 * rand(J, 1);
    
                [p_hat, c_hat, g_hat, loglike] = get_cg(X, Q,...
                    p_ini, c_ini, g_ini);
                
                p_vec(:,ii) = p_hat; c_vec(:,ii)=c_hat; g_vec(:,ii)=g_hat;
                
                llk_vec(ii) = loglike;
        end
        [~, ind] = max(llk_vec);
        
        p_dist(iter) = sum((p_vec(:,ind) - p_true(:,rr)).^2);
        c_dist(iter) = sum((c_vec(:,ind) - c_true(:,rr)).^2);
        g_dist(iter) = sum((g_vec(:,ind) - g_true(:,rr)).^2);

    %     p_mse = p_mse + sum((p_hat - p_true).^2);
    %     c_mse = c_mse + sum((c_hat - c_true).^2);
    %     g_mse = g_mse + sum((g_hat - g_true).^2);

        fprintf('Finished iteration %d,\t with MSE %1.4f, %1.4f, %1.4f \n', ...
            iter, p_dist(iter)/2^K, c_dist(iter)/J, g_dist(iter)/J);
    end
    
    p_mse(rr) = sum(p_dist)/Niter/2^K;
    c_mse(rr) = sum(c_dist)/Niter/J;
    g_mse(rr) = sum(g_dist)/Niter/J;
    fprintf('Just finished parameter set %d,\t with MSE %1.4f, %1.4f, %1.4f \n\n',...
        rr, p_mse(rr), c_mse(rr), g_mse(rr));
end


% new_name = strcat('MSE_N', num2str(N), '_Niter', num2str(Num), '.txt');
% new_file = fopen(new_name, 'w');
% fprintf(new_file, 'MSE of p: %1.10f\n', p_mse);
% fprintf(new_file, 'MSE of c: %1.10f\n', c_mse);
% fprintf(new_file, 'MSE of g: %1.10f\n', g_mse);
% fclose(new_file);


end



