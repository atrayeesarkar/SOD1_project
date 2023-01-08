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
 figure
% errorbar(x(1:n_bins),F,ebar,'.')

plot(x(1:n_bins),f,'.')


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
plot(x(1:n_bins),f,'.')
 
 
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
plot(x(1:n_bins),f,'.')
 


