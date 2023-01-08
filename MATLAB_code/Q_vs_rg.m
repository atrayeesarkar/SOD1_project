q_u = [];
for i =1:3
Q = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_',num2str(i),'/traj_comp_pbcmolcenter_WT_u.xtc.CA.Q');
% q(:,i) = load(Q);
q_u = [q_u;load(Q)];
end
rg_u = [];

for i =1:3
r = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_',num2str(i),'/rg_WT_u.txt');
rg_u = [rg_u;load(r)];
end

plot(q_u/360,rg_u(:,2),'b--','LineWidth', 2.0)
