tic
 L = 140.0;%in angstrom;
 cutoff = 15;%in angstrom;
 trj = cell(1,110);
 contacts=360;% total no. of native contacts

PUu = load('Pij_WT_uncharged_157_15_intra_single_unfolded_v3.txt');
PCu = load('Pij_WT_charged_157_15_intra_single_unfolded_v3.txt');
p1 = PCu-PUu;
[i,j]=find(triu(reshape(p1,110,110))>0.2);
k = [i,j];
[N,~] = size(k);
a = [21,23,24,30,36,40,43];
b = [53    60    61    62    66  70 71];
AA = [a,b];
d = nchoosek(AA,2);
d = [d;d(:,2),d(:,1)];
d = sortrows(d,1);
% [A,B] = meshgrid(a,b);
% c=cat(2,A',B');
% d=reshape(c,[],2);
% % % C = intersect(k,d,'rows');
% % % k = C;
k = d;
% % % % % % % % % N = N*2;
parfor i = 1:110
    t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_0/157_1/J_textfiles/','equil_',num2str(i),'.txt');
    trj{i}= load(t);
    trj{i}(1,:)=[];
end

toc
q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_0/157_1/traj_comp_pbcmolcenter_WT_c.xtc.CA.Q");
[s1,~]=size(trj{1});
[s2,~]=size(q1);
s3=s2-s1+1;
q=q1(s3:s2);
q_normalized = q/contacts;
 
unfolded_time_frame1 = find(q_normalized<0.42);
folded_time_frame1 = find(q_normalized>0.42);
parfor i = 1:110
        trj_unfolded1{i}=trj{i}(unfolded_time_frame1,:);
        trj_folded1{i}=trj{i}(folded_time_frame1,:);
end


dist = [];
parfor id1= 1:length(k)
        id1
        delta = abs(trj_unfolded1{k(id1,1)}(:,3:5)-trj_unfolded1{k(id1,2)}(:,3:5));
        delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=[dist,(sum(delta.^2,2)).^0.5];
end 
%  rg1=sqrt(sum(dist.^2,2)/(2*N));
%  figure

% Q = sum((dist<cutoff),2)/N;
%  plot(Q,rg1,'k.', 'LineWidth', 2.0)
[p,q] = find(Q(:,1)>0.79);
Q1 = Q(p,:);
dist1 = dist(p,:);
[pp,qq] = find(Q1<0.81);
Q2 = Q1(pp,:);
dist2 = dist1(pp,:);
M = mean(dist2,1);