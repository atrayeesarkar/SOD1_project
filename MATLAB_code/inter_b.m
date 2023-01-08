tic
%  L = 146.0;%in angstrom;
%  cutoff = 15;%in angstrom;
%  trj = cell(1,110);
%  contacts=360;%
% 
% PUu = load('Pij_WT_uncharged_157_15_intra_single_unfolded_v3.txt');
% PCu = load('Pij_WT_charged_157_15_intra_single_unfolded_v3.txt');
% p1 = PCu-PUu;
% [i,j]=find(triu(reshape(p1,110,110))>0.2);
% k = [i,j];
% 
%  [s,~]=size(k);
% 
% parfor i = 1:110
%     t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_20/161_1/J_textfiles/','equil_',num2str(i),'.txt');
%     trj{i}= load(t);
%    trj{i}(end-63:end,:)=[];
% end
% 
% toc
%  q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_20/161_1/Q.dat");
% % q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_20/161_1/traj_comp_pbcmolcenter.xtc.CA.Q");
% s3 = 10001;
% q=q1(s3:end,:);
% q_normalized = q/contacts;
% 
% 
% [i,j] = find(q_normalized<0.42);
% kk = [i,j];
% kk = sortrows(kk,1);
% n_unfolded = sum(q_normalized<0.42,2);
% index1 = (kk(:,1)-1)*64 +kk(:,2);
% [o,p] = find(q_normalized>0.42);
% r = [o,p];
% r = sortrows(r,1);
% n_folded = sum(q_normalized>0.42,2);
% index2 = (r(:,1)-1)*64 +r(:,2);
% parfor l= 1:110
% trj_unfolded{l} = trj{l}(index1,:);
% trj_folded{l} = trj{l}(index2,:);
% end
% 
% 
% disp("file loading done")
% dist = [];
% 
for id1=1:length(k)
    id1
            delta = abs(trj{k(id1,1)}(:,3:5)-trj{k(id1,2)}(:,3:5));
            delta(delta>L-delta)=L-delta(delta>L-delta);
            dist=[dist,(sum(delta.^2,2)).^0.5];
            
end 
 [e,f] = find(q_normalized>=0.001 & q_normalized<=0.005 );
[s,~]=size(k);
Q = sum((dist<cutoff),2)/s;
 Q = reshape(Q,64,[]);
 Q1 = Q';
 nbins = 50;
 binedges = linspace(0,1,50);
 idx = sub2ind(size(Q1),e,f); % linear indexing is the only answer 

%   h = histogram(Q1(idx),binedges, 'Normalization', 'probability');
%     plot(binedges(2:end),-log(h.Values),':b','LineWidth',4)
