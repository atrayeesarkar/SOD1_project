contacts = 360;
q = [];
n_bins = 32;
bin_width = contacts/n_bins;
bin_list = linspace(0,contacts,(contacts/bin_width)+1);
q_u=[];
for i =1:3
Q = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_',num2str(i),'/traj_comp_pbcmolcenter_WT_u.xtc.CA.Q');
% q(:,i) = load(Q);
q_u = [q_u;load(Q)];
end

% for i = 1:3
%     H = histogram(q(:,i),bin_list);
%     f(:,i) = -log(H.Values);
% end 
H = histogram(q_u,bin_list,'Normalization','probability');
f_u = -log(H.Values);
% F = mean(f,2);
x = H.BinEdges + bin_width/2;
x = x/contacts;
% j_mean = mean(jackknife(@std,f'),1);
% ebar = sqrt(sum((jackknife(@std,f')-j_mean).^2,1));
figure
subplot(2,1,1);
hold on 

axis 'square';
% ylim([]);
xlabel('q','FontSize', 30);
ylabel('\beta F(q)','FontSize', 30);
% errorbar(x(1:n_bins),F,ebar,'.')
plot(x(1:n_bins),f_u,'b--','LineWidth', 2.0)



q_m=[];
for i =1:3
Q = append('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/157_',num2str(i),'/traj_comp_pbcmolcenter_WT_m.xtc.CA.Q');
q_m = [q_m;load(Q)];
end

H = histogram(q_m,bin_list, 'Normalization', 'probability');
f_m = -log(H.Values);
% F = mean(f,2);
x = H.BinEdges + bin_width/2;
x = x/contacts;
% j_mean = mean(jackknife(@std,f'),1);
% ebar = sqrt(sum((jackknife(@std,f')-j_mean).^2,1));
% errorbar(x(1:n_bins),F,ebar,'.-')
plot(x(1:n_bins),f_m,'r.-', 'LineWidth', 2.0)




q_c = [];
for i =1:3
Q = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/157_',num2str(i),'/traj_comp_pbcmolcenter_WT_c.xtc.CA.Q');
q_c = [q_c;load(Q)];
end

% for i = 1:3
%     H = histogram(q(:,i),bin_list);
%     f(:,i) = -log(H.Values);
% end 
H = histogram(q_c,bin_list, 'Normalization', 'probability');
f_c = -log(H.Values);
% F = mean(f,2);
x = H.BinEdges + bin_width/2;
x = x/contacts;
% j_mean = mean(jackknife(@std,f'),1);
% ebar = sqrt(sum((jackknife(@std,f')-j_mean).^2,1));
% errorbar(x(1:n_bins),F,ebar,'*')
plot(x(1:n_bins),f_c,'k:','LineWidth', 2.0)

legend('WT_u','Mutant','WT_c', 'FontSize', 30 );

%%%%%%%%%%%%%%% crowded cases %%%%%%%%%%%%%%%

contacts = 360;
q = [];
n_bins = 32;
bin_width = contacts/n_bins;
bin_list = linspace(0,contacts,(contacts/bin_width)+1);
n = 64;
% f=[];
% parfor i =1:n
% Q = append('/Users/atrayee/Desktop/SOD_project/Data/WT_uncharged/phic_20/161/Q/','q_',num2str(i),'.dat');
% q(:,i) = load(Q);
% H = histogram(q(:,i),bin_list);
% h =[h; H.Values];
% f(:,i) = -log(H.Values);
% end
Q = load('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_20/161_1/traj_comp_pbcmolcenter_WT_u.xtc.CA.Q');
H = histogram(Q,bin_list,'Normalization', 'probability');
f = -log(H.Values);
% [p,q]=find((isinf(f)));
% f(:,q) = [];
% H = histogram(q(:,1),bin_list);


% F = -log(mean(h,1));
x = H.BinEdges + bin_width/2;
x = x/contacts;
% j_mean = mean(jackknife(@std,f'),1);
% ebar = sqrt(sum((jackknife(@std,f')-j_mean).^2,1));
 subplot(2,1,2);
 hold on 

axis 'square';
% ylim([]);
xlabel('q','FontSize', 30);
ylabel('\beta F(q)','FontSize', 30);

% errorbar(x(1:n_bins),F,ebar,'.')

plot(x(1:n_bins),f,'b--','LineWidth', 2.0)


 hold on 
 
 % f=[];
% q = [];
% parfor i =1:n
% Q = append('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_20/161/Q/','q_',num2str(i),'.dat');
% q(:,i) = load(Q);
% H = histogram(q(:,i),bin_list,'Normalization', 'probability');
% h =[h; H.Values];
% f(:,i) = -log(H.Values);
% end
% % [p,q]=find((isinf(f)));
% % f(:,q) = [];
% 
% F = -log(mean(h,1));
% x = H.BinEdges + bin_width/2;
% x = x/contacts;
% j_mean = mean(jackknife(@std,f'),1);
% ebar = sqrt(sum((jackknife(@std,f')-j_mean).^2,1));
% errorbar(x(1:n_bins),F,ebar,'.-')
Q = load('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_20/161_1/traj_comp_pbcmolcenter_WT_m.xtc.CA.Q');
H = histogram(Q,bin_list,'Normalization', 'probability');
f = -log(H.Values);
plot(x(1:n_bins),f,'r-.','LineWidth', 2.0)
 
 
% q = [];
% f=[];
% h = [];
% parfor i =1:n
% Q = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_20/161/Q/','q_',num2str(i),'.dat');
% q(:,i) = load(Q);
% H = histogram(q(:,i),bin_list);
% h =[h; H.Values];
%  f(:,i) = -log(H.Values);
% end
% % [p,q]=find((isinf(f)));
% % f(:,q) = [];
% H = histogram(q(:,1),bin_list); 
Q = load('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_20/161_1/traj_comp_pbcmolcenter_WT_c.xtc.CA.Q');
H = histogram(Q,bin_list,'Normalization', 'probability');
f = -log(H.Values);
% F = -log(mean(h,1));
% x = H.BinEdges + bin_width/2;
% x = x/contacts;
% j_mean = mean(jackknife(@std,f'),1);
% ebar = sqrt(sum((jackknife(@std,f')-j_mean).^2,1));
% errorbar(x(1:n_bins),F,ebar,'*')
plot(x(1:n_bins),f,'k:','LineWidth', 2.0)
 
legend('WT_u','Mutant','WT_c', 'FontSize', 30 );


