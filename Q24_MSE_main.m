%%%%%%%%%%%%%%%%% Q24: random generation and generic identifiability %%%%%%%%%%%%%%%%%%%
Q = [1 0; 0 1; 1 0; 0 1];
[J, K] = size(Q);

I_use = get_I(2, J);
I_full = get_I(J,J);
I_full1 = [zeros(1,J); I_full];

A = binary(0:(2^K-1), K);

% randomly generate Nset sets of paraemters
Nset = 100;
p_true = zeros(2^K, Nset);
c_true = zeros(J, Nset);
g_true = zeros(J, Nset);

rng('default');
for rr = 1:Nset
    % pn = 0.1 * rand(4,1) + 0.1;
    pn = gamrnd(5*ones(1,2^K), ones(1,2^K));
    p_true(:,rr) = pn / sum(pn);
    c_true(:,rr) = 0.9 - 0.2 * rand(J, 1);
    g_true(:,rr) = 0.1 + 0.2 * rand(J, 1);
end



[p_mse50, c_mse50, g_mse50] = get_MSE(50,p_true, c_true, g_true);



p_mse_mat = [p_mse100, p_mse1e3, p_mse1e4, p_mse1e5];

c_mse_mat = [c_mse100, c_mse1e3, c_mse1e4, c_mse1e5];

g_mse_mat = [g_mse100, g_mse1e3, g_mse1e4, g_mse1e5];


figure
h = boxplot(p_mse_mat, 'Labels', {'1e2', '1e3', '1e4', '1e5'});
set(h,{'linew'},{2})
grid on
hold on
plot(1:4, median(p_mse_mat), '-x', 'LineWidth', 3, 'MarkerSize', 12);
xlabel('N')
pbaspect([1 1 1]);
set(gca,'fontsize',15)
print('-r400', 'p_mse1', '-dpng');



fig = figure;
h = boxplot(c_mse_mat, 'Labels', {'1e2', '1e3', '1e4', '1e5'});
set(h,{'linew'},{2})
xlabel('N')
grid on
hold on
plot(1:4, median(c_mse_mat), '-x', 'LineWidth', 3, 'MarkerSize', 12);
pbaspect([1 1 1]);
set(gca,'fontsize',15)
print('-r400', 's_mse1', '-dpng');



figure
h = boxplot(g_mse_mat, 'Labels',  {'1e2', '1e3', '1e4', '1e5'});
set(h,{'linew'},{2})
xlabel('N')
grid on
hold on
plot(1:4, median(g_mse_mat), '-x', 'LineWidth', 3, 'MarkerSize', 12);
pbaspect([1 1 1]);
set(gca,'fontsize',15)
print('-r400', 'g_mse1', '-dpng');



median(p_mse_mat)
median(c_mse_mat)
median(g_mse_mat)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MSE_100.mat')
load('MSE_1e3.mat')
load('MSE_1e4.mat')
load('MSE_1e5.mat')

c_mse1e5;

% hOutliers = findobj(bp_x,'Tag','Outliers');
% yy = get(hOutliers,'YData');
% ind_ot_c = find(c_mse1e5 == yy);

ind_ot_c = find(c_mse1e5>0.005);
p_ot = p_true(:,ind_ot_c);

[p_ot(2,:)./p_ot(1,:); p_ot(4,:)./p_ot(3,:)]

[p_ot(3,:)./p_ot(1,:); p_ot(4,:)./p_ot(2,:)]



g_mse1e5;
ind_ot_g = find(g_mse1e5>0.005)
p_ot = p_true(:,ind_ot_g);

[p_ot(2,:)./p_ot(1,:); p_ot(4,:)./p_ot(3,:)]
figure
plot(p_ot(2,:)./p_ot(1,:), p_ot(4,:)./p_ot(3,:),'*')

[p_ot(3,:)./p_ot(1,:); p_ot(4,:)./p_ot(2,:)]



p_mse1e5;

p_quant = quantile(p_mse1e5, [0 0.25 0.5 0.75 1]);
thres = p_quant(4) + (p_quant(4)-p_quant(2))*1.57/100;

ind_ot_p = find(p_mse1e5>thres);
int_reg = find(p_mse1e5 <= thres);

p_ot = p_true(:, ind_ot_p);
p_reg = p_true(:,p_mse1e5 <= thres);


[p_ot(2,:)./p_ot(1,:); p_ot(4,:)./p_ot(3,:)]
[p_ot(3,:)./p_ot(1,:); p_ot(4,:)./p_ot(2,:)]

[p_true(2,:)./p_true(1,:); p_true(4,:)./p_true(3,:)]
[p_true(3,:)./p_true(1,:); p_true(4,:)./p_true(2,:)]



% figure
% plot(p_ot(2,:)./p_ot(1,:), p_ot(4,:)./p_ot(3,:),'*')
% hold on
% xlim([0, 1]); ylim([0, 1])
% refline(1,0)
% pbaspect([1 1 1]);
% xlabel('p_{01}/p_{00}'); ylabel('p_{11}/p_{10}')
% plot(p_reg(2,:)./p_reg(1,:), p_reg(4,:)./p_reg(3,:), '+')
% 
% figure
% plot(p_reg(3,:)./p_reg(1,:), p_reg(4,:)./p_reg(2,:), '+')
% hold on
% plot(p_ot(3,:)./p_ot(1,:), p_ot(4,:)./p_ot(2,:),'*')
% xlim([0, 1]); ylim([0, 1])
% refline(1,0)
% pbaspect([1 1 1]);
% xlabel('p_{10}/p_{00}'); ylabel('p_{11}/p_{01}')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(p_reg(2,:)./p_reg(1,:), p_reg(4,:)./p_reg(3,:), '+')
hold on
plot(p_ot(2,:)./p_ot(1,:), p_ot(4,:)./p_ot(3,:),'*')
xlim([0, 10]); ylim([0, 10])
refline(1,0)
pbaspect([1 1 1]);
xlabel('p_{01}/p_{00}'); ylabel('p_{11}/p_{10}')
set(gca,'fontsize',15)
print('-r400', 'ratio1', '-dpng');



figure
plot(p_reg(3,:)./p_reg(1,:), p_reg(4,:)./p_reg(2,:), '+')
hold on
plot(p_ot(3,:)./p_ot(1,:), p_ot(4,:)./p_ot(2,:),'*')
xlim([0, 10]); ylim([0, 10])
refline(1,0)
pbaspect([1 1 1]);
xlabel('p_{10}/p_{00}'); ylabel('p_{11}/p_{01}')
set(gca,'fontsize',15)
print('-r400', 'ratio2', '-dpng');





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% prod_reg = (p_reg(2,:)./p_reg(1,:) -  p_reg(4,:)./p_reg(3,:))...
%          .* (p_reg(3,:)./p_reg(1,:) -  p_reg(4,:)./p_reg(2,:));
%          
% prod_ot = (p_ot(2,:)./p_ot(1,:) -  p_ot(4,:)./p_ot(3,:))...
%          .* (p_ot(3,:)./p_ot(1,:) -  p_ot(4,:)./p_ot(2,:));
%      
% quantile(abs(prod_reg), [0 0.25 0.5 0.75])
% quantile(abs(prod_ot), [0 0.25 0.5 0.75])
% 
% boxplot(prod_reg)
%      
% figure
% plot(ind_reg(abs(prod_reg)<40), prod_reg(abs(prod_reg)<40), '+')
% hold on
% plot(ind_ot_p, prod_ot, '*')
% 
% 
% figure
% plot(p_reg(2,:)./p_reg(1,:), p_reg(4,:)./p_reg(3,:), '+')
% hold on
% plot(p_ot(2,:)./p_ot(1,:), p_ot(4,:)./p_ot(3,:),'*')
% xlim([0, 10]); ylim([0, 10])
% refline(1,0)
% set(gca,'fontsize',15)
% 
% hold on
% 
% plot(p_reg(3,:)./p_reg(1,:), p_reg(4,:)./p_reg(2,:), '+')
% hold on
% plot(p_ot(3,:)./p_ot(1,:), p_ot(4,:)./p_ot(2,:),'*')
% xlim([0, 10]); ylim([0, 10])
% refline(1,0)
% pbaspect([1 1 1]);
% 


p_quant = quantile(p_mse1e5, [0 0.25 0.5 0.75 1]);
thres = p_quant(4) + (p_quant(4)-p_quant(2))*1.57/100;

[ b, ix ] = sort( p_mse1e5, 'descend' );

ind_ot = ix(1:20);
ind_reg = find(p_mse1e5 < p_mse1e5(ix(20)));

p_ot = p_true(:, ind_ot);
p_reg = p_true(:, ind_reg);


reg_ratio = (p_reg(1,:).*p_reg(4,:)) ./ (p_reg(2,:).*p_reg(3,:));
ot_ratio = (p_ot(1,:).*p_ot(4,:)) ./ (p_ot(2,:).*p_ot(3,:));
figure
plot( (p_reg(1,:).*p_reg(4,:)) ,  (p_reg(2,:).*p_reg(3,:)) , '+', 'MarkerSize', 15)
hold on
plot( (p_ot(1,:).*p_ot(4,:)) ,  (p_ot(2,:).*p_ot(3,:)) , '*', 'MarkerSize', 15)
%xlim([0, 10]); ylim([0, 10])
refline(1,0)
pbaspect([1 1 1]);
xlabel('p_{00}p_{11}'); ylabel('p_{01}p_{10}')
set(gca,'fontsize',20)

% %%%% new color map %%%%
% reg_ratio = (p_true(1,:).*p_true(4,:)) ./ (p_true(2,:).*p_true(3,:));
% ot_ratio = (p_ot(1,:).*p_ot(4,:)) ./ (p_ot(2,:).*p_ot(3,:));
% figure
% % plot( (p_true(1,:).*p_true(4,:)) ,  (p_true(2,:).*p_true(3,:)) , '+', 'MarkerSize', 10)
% % hold on
% % plot( (p_ot(1,:).*p_ot(4,:)) ,  (p_ot(2,:).*p_ot(3,:)) , '*', 'MarkerSize', 10)
% refline(1,0)
% pbaspect([1 1 1]);
% xlabel('p_{00}p_{11}'); ylabel('p_{01}p_{10}')
% set(gca,'fontsize',20)
% x = p_true(1,:).*p_true(4,:);
% y = (p_true(2,:).*p_true(3,:));
% c = p_mse1e5;
% scatter(x,y,[],c)
% colormap(hot(64));
% colorbar;
% %%%% new color map %%%%


print('-r400', 'ratio20', '-dpng');




%%% try %%%
% randomly generate Nset sets of parameters
J = 4; K = 2;
Nset = 100;
p_true = zeros(2^K, Nset);
c_true = zeros(J, Nset);
g_true = zeros(J, Nset);

rng('default');
for rr = 1:Nset
    % pn = 0.1 * rand(4,1) + 0.1;
    pn = gamrnd(3*ones(1,2^K), ones(1,2^K));
    p_true(:,rr) = pn / sum(pn);
    c_true(:,rr) = 0.9 - 0.2 * rand(J, 1);
    g_true(:,rr) = 0.1 + 0.2 * rand(J, 1);
end

mean(min(p_true))
sqrt(var(min(p_true)))
min(min(p_true))

figure
plot( (p_true(1,:).*p_true(4,:)), (p_true(2,:).*p_true(3,:)), '+', 'MarkerSize', 15)
hold on
refline(1,0)
pbaspect([1 1 1]);
xlabel('p_{00}\times p_{11}'); ylabel('p_{01}\times p_{10}')
%xlim([0 0.25]); ylim([0 0.25])
set(gca,'fontsize',20)

