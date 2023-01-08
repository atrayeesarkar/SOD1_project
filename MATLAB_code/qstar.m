% tic
%  L = 140.0;%in angstrom;
%  cutoff = 15;%in angstrom;
%  trj = cell(1,110);
%  contacts=360;% total no. of native contacts
% 
% PUu = load('Pij_WT_uncharged_157_15_intra_single_unfolded_v2.txt');
% PCu = load('Pij_WT_charged_157_15_intra_single_unfolded_v2.txt');
% p1 = PCu-PUu;
% [i,j]=find(triu(reshape(p1,110,110))>0.2);
% k = [i,j];
% c = load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Uncharged/phic_0/157_1/con.txt");
% c = c(:,1:2);
% native = intersect(k,c,'rows');
% non_native = setdiff(k,native,'rows');
% [s,~]=size(k);
% 
% parfor i = 1:110
%     t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_0/157_1/J_textfiles/','equil_',num2str(i),'.txt');
%     trj{i}= load(t);
%     trj{i}(1,:)=[];
% end
% 
% toc
% q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_0/157_1/traj_comp_pbcmolcenter_WT_m.xtc.CA.Q");
% [s1,~]=size(trj{1});
% [s2,~]=size(q1);
% s3=s2-s1+1;
% q=q1(s3:s2);
% q_normalized = q/contacts;
%  
% unfolded_time_frame1 = find(q_normalized<0.42);
% folded_time_frame1 = find(q_normalized>0.42);
% parfor i = 1:110
%         trj_unfolded1{i}=trj{i}(unfolded_time_frame1,:);
%         trj_folded1{i}=trj{i}(folded_time_frame1,:);
% end
% 
% 
% dist = [];
% parfor id1= 1:length(k)
%         id1
%         delta = abs(trj{k(id1,1)}(:,3:5)-trj{k(id1,2)}(:,3:5));
%         delta(delta>L-delta)=L-delta(delta>L-delta);
%         dist=[dist,(sum(delta.^2,2)).^0.5];
% end 
% 
%  
%  Q = sum((dist<cutoff),2)/s;
%  [e,f] = find(q_normalized>=0.14 & q_normalized<=0.16 );
%   figure
%    nbins = 50;
%    binedges = linspace(0,1,50);
%   h = histogram(Q(e),binedges,'Normalization', 'probability','LineWidth',4);
%  h.FaceColor = 'None';
%  h.EdgeColor = 'r';
%  h.LineStyle = '-.';
%  h.LineWidth = 2;
%  h.FaceColor = 'None';
%  nbins = 50;
%  binedges = linspace(0,1,50);
%  figure
%  h = histogram(Q,binedges);
%  
%    plot(binedges(2:end),-log(h.Values),':b','LineWidth',4)
%  figure
%  h = histogram2(q_normalized,Q,'XBinEdges',binedges, 'YBinEDges', binedges,'Normalization', 'probability');
%  s=surf(binedges(2:end),binedges(2:end),-log(h.Values)');
%  s.EdgeColor = 'None';
 
%  dist = [];
% parfor id1= 1:length(k)
%         id1
%         delta = abs(trj_unfolded1{k(id1,1)}(:,3:5)-trj_unfolded1{k(id1,2)}(:,3:5));
%         delta(delta>L-delta)=L-delta(delta>L-delta);
%         dist=[dist,(sum(delta.^2,2)).^0.5];
% end 
% [s,~]=size(k);
%  
%  Qu = sum((dist<cutoff),2)/s;
%  figure
%  xbinedges = linspace(0,0.42,25);
%  ybinedges = linspace(0,1,25);
%  h = histogram(Qu,ybinedges,'Normalization', 'probability');
%  h.FaceColor = 'None';
%  h.EdgeColor = 'b';
%  h.LineStyle = '--';
%  h.LineWidth = 2;
%  
%  
%  
%  
%  
% parfor i = 1:110
%     t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_0/157_1/J_textfiles/','equil_',num2str(i),'.txt');
%     trj{i}= load(t);
%     trj{i}(1,:)=[];
% end
% 
% toc
% q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_0/157_1/traj_comp_pbcmolcenter_WT_m.xtc.CA.Q");
% [s1,~]=size(trj{1});
% [s2,~]=size(q1);
% s3=s2-s1+1;
% q=q1(s3:s2);
% q_normalized = q/contacts;
%  
% unfolded_time_frame1 = find(q_normalized<0.42);
% folded_time_frame1 = find(q_normalized>0.42);
% parfor i = 1:110
%         trj_unfolded1{i}=trj{i}(unfolded_time_frame1,:);
%         trj_folded1{i}=trj{i}(folded_time_frame1,:);
% end
% 
% 
% dist = [];
% parfor id1= 1:length(k)
%         id1
%         delta = abs(trj{k(id1,1)}(:,3:5)-trj{k(id1,2)}(:,3:5));
%         delta(delta>L-delta)=L-delta(delta>L-delta);
%         dist=[dist,(sum(delta.^2,2)).^0.5];
% end 
% 
%  
%  Q = sum((dist<cutoff),2)/s;
% %  figure
% %  nbins = 25;
% %  binedges = linspace(0,1,25);
% %  h = histogram(Q,binedges,'Normalization', 'probability');
% %  h.FaceColor = 'None';
% %  h.EdgeColor = 'r';
% %  h.LineStyle = '-.';
% %  h.LineWidth = 2;
% %  h.FaceColor = 'None';
% %  nbins = 50;
% %  binedges = linspace(0,1,50);
% %  figure
% %  h = histogram(Q,binedges);
%  
% %  plot(binedges(2:end),-log(h.Values))
%  figure
%  h = histogram2(q_normalized,Q,'XBinEdges',binedges, 'YBinEDges', binedges,'Normalization', 'probability');
%  s=surf(binedges(2:end),binedges(2:end),-log(h.Values)');
%  s.EdgeColor = 'None';
%  
%  dist = [];
% parfor id1= 1:length(k)
%         id1
%         delta = abs(trj_unfolded1{k(id1,1)}(:,3:5)-trj_unfolded1{k(id1,2)}(:,3:5));
%         delta(delta>L-delta)=L-delta(delta>L-delta);
%         dist=[dist,(sum(delta.^2,2)).^0.5];
% end 
% [s,~]=size(k);
%  
%  Qu = sum((dist<cutoff),2)/s;
%  figure
%  xbinedges = linspace(0,0.42,25);
%  ybinedges = linspace(0,1,25);
%  h = histogram(Qu,ybinedges,'Normalization', 'probability');
%  h.FaceColor = 'None';
%  h.EdgeColor = 'r';
%  h.LineStyle = '-.';
%  h.LineWidth = 2;
% 
% %  xbinedges = linspace(0,0.42,50);
% %  ybinedges = linspace(0,1,50);
% %  figure
% %  h = histogram(Qu,ybinedges);
% %  plot(ybinedges(2:end),-log(h.Values))
% %  figure
% %  h = histogram2(q_normalized(unfolded_time_frame1),Qu,'XBinEdges',xbinedges, 'YBinEDges', ybinedges,'Normalization', 'probability');
% %  s = surf(xbinedges(2:end),ybinedges(2:end),-log(h.Values)');
% %  s.EdgeColor = 'None';
%  
% %  dist = [];
% % parfor id1= 1:length(k)
% %         id1
% %         delta = abs(trj_folded1{k(id1,1)}(:,3:5)-trj_folded1{k(id1,2)}(:,3:5));
% %         delta(delta>L-delta)=L-delta(delta>L-delta);
% %         dist=[dist,(sum(delta.^2,2)).^0.5];
% % end 
% % 
% %  [s,~]=size(k);
% %  Qf = sum((dist<cutoff),2)/s;
% %  figure
% %  xbinedges = linspace(0.42,1,50);
% %  ybinedges = linspace(0,1,50);
% %  h = histogram(Qf,ybinedges,'Normalization', 'probability');
% %  h.FaceColor = 'None';
% %  figure
% %  h = histogram(Qf,ybinedges);
% %  plot(ybinedges(2:end),-log(h.Values))
% %  figure
% %  h = histogram2(q_normalized(folded_time_frame1),Qf,'XBinEdges',xbinedges, 'YBinEDges', ybinedges,'Normalization', 'probability');
% %  s = surf(xbinedges(2:end),ybinedges(2:end),-log(h.Values)');
% %  s.EdgeColor = 'None';
% 
% 
% parfor i = 1:110
%     t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_0/157_1/J_textfiles/','equil_',num2str(i),'.txt');
%     trj{i}= load(t);
%     trj{i}(1,:)=[];
% end
% 
% toc
% q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_0/157_1/traj_comp_pbcmolcenter_WT_c.xtc.CA.Q");
% [s1,~]=size(trj{1});
% [s2,~]=size(q1);
% s3=s2-s1+1;
% q=q1(s3:s2);
% q_normalized = q/contacts;
%  
% unfolded_time_frame1 = find(q_normalized<0.42);
% folded_time_frame1 = find(q_normalized>0.42);
% parfor i = 1:110
%         trj_unfolded1{i}=trj{i}(unfolded_time_frame1,:);
%         trj_folded1{i}=trj{i}(folded_time_frame1,:);
% end
% 
% 
% dist = [];
% parfor id1= 1:length(k)
%         id1
%         delta = abs(trj{k(id1,1)}(:,3:5)-trj{k(id1,2)}(:,3:5));
%         delta(delta>L-delta)=L-delta(delta>L-delta);
%         dist=[dist,(sum(delta.^2,2)).^0.5];
% end 
% 
%  
%  Q = sum((dist<cutoff),2)/s;
% %  figure
% %  nbins = 25;
% %  binedges = linspace(0,1,25);
% %  h = histogram(Q,binedges,'Normalization', 'probability');
% %  h.FaceColor = 'None';
% %  h.EdgeColor = 'r';
% %  h.LineStyle = '-.';
% %  h.LineWidth = 2;
% %  h.FaceColor = 'None';
% %  nbins = 50;
% %  binedges = linspace(0,1,50);
% %  figure
% %  h = histogram(Q,binedges);
%  
% %  plot(binedges(2:end),-log(h.Values))
%  figure
%  h = histogram2(q_normalized,Q,'XBinEdges',binedges, 'YBinEDges', binedges,'Normalization', 'probability');
%  s=surf(binedges(2:end),binedges(2:end),-log(h.Values)');
%  s.EdgeColor = 'None';
%  
%  dist = [];
% parfor id1= 1:length(k)
%         id1
%         delta = abs(trj_unfolded1{k(id1,1)}(:,3:5)-trj_unfolded1{k(id1,2)}(:,3:5));
%         delta(delta>L-delta)=L-delta(delta>L-delta);
%         dist=[dist,(sum(delta.^2,2)).^0.5];
% end 
% [s,~]=size(k);
%  
%  Qu = sum((dist<cutoff),2)/s;
%  figure
%  xbinedges = linspace(0,0.42,25);
%  ybinedges = linspace(0,1,25);
%  h = histogram(Qu,ybinedges,'Normalization', 'probability');
%  h.FaceColor = 'None';
%  h.EdgeColor = 'k';
%  h.LineStyle = ':';
%  h.LineWidth = 2;
% 
% 
%  toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






tic
 L = 146.0;%in angstrom;
 cutoff = 15;%in angstrom;
 trj = cell(1,110);
 contacts=360;%

PUu = load('Pij_WT_uncharged_157_15_intra_single_unfolded_v3.txt');
PCu = load('Pij_WT_charged_157_15_intra_single_unfolded_v3.txt');
p1 = PCu-PUu;
[i,j]=find(triu(reshape(p1,110,110))>0.2);
k = [i,j];

 [s,~]=size(k);

parfor i = 1:110
    t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_20/161_1/J_textfiles/','equil_',num2str(i),'.txt');
    trj{i}= load(t);
   trj{i}(end-63:end,:)=[];
end

toc
 q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_20/161_1/Q.dat");
% q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_20/161_1/traj_comp_pbcmolcenter.xtc.CA.Q");
s3 = 10001;
q=q1(s3:end,:);
q_normalized = q/contacts;


[i,j] = find(q_normalized<0.42);
kk = [i,j];
kk = sortrows(kk,1);
n_unfolded = sum(q_normalized<0.42,2);
index1 = (kk(:,1)-1)*64 +kk(:,2);
[o,p] = find(q_normalized>0.42);
r = [o,p];
r = sortrows(r,1);
n_folded = sum(q_normalized>0.42,2);
index2 = (r(:,1)-1)*64 +r(:,2);
parfor l= 1:110
trj_unfolded{l} = trj{l}(index1,:);
trj_folded{l} = trj{l}(index2,:);
end


disp("file loading done")
dist = [];
% 
parfor id1=1:length(k)
    id1
            delta = abs(trj{k(id1,1)}(:,3:5)-trj{k(id1,2)}(:,3:5));
            delta(delta>L-delta)=L-delta(delta>L-delta);
            dist=[dist,(sum(delta.^2,2)).^0.5];
            
end 
 [e,f] = find(q_normalized>=0.14 & q_normalized<=0.16 );
% H = zeros(64,49);
% H1 = zeros(64,49,49);
[s,~]=size(k);
Q = sum((dist<cutoff),2)/s;
 Q = reshape(Q,64,[]);
 Q1 = Q';
 nbins = 50;
 binedges = linspace(0,1,50);
 idx = sub2ind(size(Q1),e,f); % linear indexing is the only answer 
figure
% l = 1;
% for i = 1:64
%  
  h = histogram(Q1(idx),binedges, 'Normalization', 'probability');
    plot(binedges(2:end),-log(h.Values),':b','LineWidth',4)
%  H(i,:) = h.Values; 
% figure
% q_normalized=reshape(q_normalized',[],1);
%   h = histogram2(q_normalized,Q1,'XBinEdges',binedges, 'YBinEDges', binedges,'Normalization', 'probability');
%  H1(i,:,:) = h1.Values; 
%  l = l+64;
% end
% P = mean(H,1);
% P1 = mean(H1,1);
% P3(:,:) = P1(1,:,:);
% plot(binedges(2:end),P)
% figure
%   s=surf(binedges(2:end),binedges(2:end),-log(h.Values)');
%   s.EdgeColor = 'None';
%   
 

% dist = [];
% 
% parfor id1=1:length(k)
%     id1
%             delta = abs(trj_unfolded{k(id1,1)}(:,3:5)-trj_unfolded{k(id1,2)}(:,3:5));
%             delta(delta>L-delta)=L-delta(delta>L-delta);
%             dist=[dist,(sum(delta.^2,2)).^0.5];
%             
% end 
% % H = zeros(max(n_unfolded),49);
% [s,~]=size(k);
% Qu = sum((dist<cutoff),2)/s;
% ybinedges = linspace(0,1,25);
% figure
% % l = 1;
% % nn = nan(max(n_unfolded),10001);
% % for i = 1:length(n_unfolded)
% %     nn(1:n_unfolded(i),i) = Qu(l:l+n_unfolded(i)-1);
% %      l = l+n_unfolded(i);
% % end 
% 
% % for i = 1:max(n_unfolded)
% %  h = histogram(nn(i,:)',ybinedges, 'Normalization', 'probability');
% %  H(i,:) = h.Values; 
% % %  l = l+n_unfolded(i);
% % % end
% % % P = mean(H,1);
% % plot(ybinedges(2:end),P)
% h = histogram(Qu,ybinedges, 'Normalization', 'probability');
% h.FaceColor = 'None';
% h.EdgeColor = 'b';
% h.LineStyle = '--';
% h.LineWidth = 2;
% 
% % dist = [];
% % 
% % parfor id1=1:length(k)
% %     id1
% %             delta = abs(trj_folded{k(id1,1)}(:,3:5)-trj_folded{k(id1,2)}(:,3:5));
% %             delta(delta>L-delta)=L-delta(delta>L-delta);
% %             dist=[dist,(sum(delta.^2,2)).^0.5];
% %             
% % end 
% % % H = zeros(64,49);
% % % [s,~]=size(k);
% % % Qf = sum((dist<cutoff),2)/s;
% % % xbinedges = linspace(0,0.42,50);
% % % ybinedges = linspace(0,1,50);
% % % figure
% % % l = 1;
% % % for i = 1:length(n_folded)
% % %  h = histogram(Qf(l:l+n_folded(i)-1)',ybinegdes, 'Normalization', 'probability');
% % %  H(i,:) = h.Values; 
% % %  l = l+n_folded(i);
% % % end
% % % P = mean(H,1);
% % % plot(ybinedges(2:end),P)
% % H = zeros(max(n_folded),49);
% % [s,~]=size(k);
% % Qu = sum((dist<cutoff),2)/s;
% % ybinedges = linspace(0,1,50);
% % figure
% % l = 1;
% % nn = nan(max(n_folded),10001);
% % for i = 1:length(n_folded)
% %     nn(1:n_folded(i),i) = Qu(l:l+n_folded(i)-1);
% %      l = l+n_folded(i);
% % end 
% % 
% % for i = 1:max(n_folded)
% %  h = histogram(nn(i,:)',ybinedges, 'Normalization', 'probability');
% %  H(i,:) = h.Values; 
% %  l = l+n_unfolded(i);
% % end
% % P = mean(H,1);
% % plot(ybinedges(2:end),P)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% parfor i = 1:110
%     t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_20/161_1/J_textfiles/','equil_',num2str(i),'.txt');
%     trj{i}= load(t);
%    trj{i}(end-63:end,:)=[];
% end
% q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/WT_Charged/phic_20/161_1/Q.dat");
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
% %  disp("file loading done")
% %  dist = [];
% % % 
% % parfor id1=1:length(k)
% %     id1
% %             delta = abs(trj{k(id1,1)}(:,3:5)-trj{k(id1,2)}(:,3:5));
% %             delta(delta>L-delta)=L-delta(delta>L-delta);
% %             dist=[dist,(sum(delta.^2,2)).^0.5];
% %             
% % end 
% % % H = zeros(64,49);
% % % H1 = zeros(64,49,49);
% % % 
% %  Q = sum((dist<cutoff),2)/s;
% % % Q = reshape(Q,64,[]);
% %  nbins = 50;
% %  binedges = linspace(0,1,50);
% %  figure
% % % l = 1;
% % % for i = 1:64
% % %  h = histogram(Q(i,:)',binedges, 'Normalization', 'probability');
% % %  H(i,:) = h.Values; 
% %   h = histogram2(q_normalized(:,i),Q','XBinEdges',binedges, 'YBinEDges', binedges,'Normalization', 'probability');
% % %  H1(i,:,:) = h1.Values; 
% % %  l = l+64;
% % % end
% % % P = mean(H,1);
% % % P1 = mean(H1,1);
% % % P3(:,:) = P1(1,:,:);
% % % plot(binedges(2:end),P)
% % % figure
% %   s=surf(binedges(2:end),binedges(2:end),-log(h.Values)');
% %   s.EdgeColor = 'None';
% % 
% 
% dist = [];
% 
% parfor id1=1:length(k)
%     id1
%             delta = abs(trj_unfolded{k(id1,1)}(:,3:5)-trj_unfolded{k(id1,2)}(:,3:5));
%             delta(delta>L-delta)=L-delta(delta>L-delta);
%             dist=[dist,(sum(delta.^2,2)).^0.5];
%             
% end 
% % H = zeros(max(n_unfolded),49);
% [s,~]=size(k);
% Qu = sum((dist<cutoff),2)/s;
% ybinedges = linspace(0,1,25);
% figure
% % l = 1;
% % nn = nan(max(n_unfolded),10001);
% % for i = 1:length(n_unfolded)
% %     nn(1:n_unfolded(i),i) = Qu(l:l+n_unfolded(i)-1);
% %      l = l+n_unfolded(i);
% % end 
% 
% % for i = 1:max(n_unfolded)
% %  h = histogram(nn(i,:)',ybinedges, 'Normalization', 'probability');
% %  H(i,:) = h.Values; 
% % %  l = l+n_unfolded(i);
% % % end
% % % P = mean(H,1);
% % plot(ybinedges(2:end),P)
% h = histogram(Qu,ybinedges, 'Normalization', 'probability');
% h.FaceColor = 'None';
% h.EdgeColor = 'k';
% h.LineStyle = ':';
% h.LineWidth = 2;
% 
% % dist = [];
% % 
% % parfor id1=1:length(k)
% %     id1
% %             delta = abs(trj_folded{k(id1,1)}(:,3:5)-trj_folded{k(id1,2)}(:,3:5));
% %             delta(delta>L-delta)=L-delta(delta>L-delta);
% %             dist=[dist,(sum(delta.^2,2)).^0.5];
% %             
% % end 
% % % H = zeros(64,49);
% % % [s,~]=size(k);
% % % Qf = sum((dist<cutoff),2)/s;
% % % xbinedges = linspace(0,0.42,50);
% % % ybinedges = linspace(0,1,50);
% % % figure
% % % l = 1;
% % % for i = 1:length(n_folded)
% % %  h = histogram(Qf(l:l+n_folded(i)-1)',ybinegdes, 'Normalization', 'probability');
% % %  H(i,:) = h.Values; 
% % %  l = l+n_folded(i);
% % % end
% % % P = mean(H,1);
% % % plot(ybinedges(2:end),P)
% % H = zeros(max(n_folded),49);
% % [s,~]=size(k);
% % Qu = sum((dist<cutoff),2)/s;
% % ybinedges = linspace(0,1,50);
% % figure
% % l = 1;
% % nn = nan(max(n_folded),10001);
% % for i = 1:length(n_folded)
% %     nn(1:n_folded(i),i) = Qu(l:l+n_folded(i)-1);
% %      l = l+n_folded(i);
% % end 
% % 
% % for i = 1:max(n_folded)
% %  h = histogram(nn(i,:)',ybinedges, 'Normalization', 'probability');
% %  H(i,:) = h.Values; 
% %  l = l+n_unfolded(i);
% % end
% % P = mean(H,1);
% % plot(ybinedges(2:end),P)
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % subplot(2,1,2)
% % hold on
% % axis 'square';
% % xlabel('b','FontSize', 30);
% % ylabel('probability','FontSize', 30);
% 
% parfor i = 1:110
%     t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_20/161_1/J_textfiles/','equil_',num2str(i),'.txt');
%     trj{i}= load(t);
%    trj{i}(end-63:end,:)=[];
% end
% q1 =load("/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_20/161_1/Q.dat");
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
% parfor id1=1:length(k)
%     id1
%             delta = abs(trj{k(id1,1)}(:,3:5)-trj{k(id1,2)}(:,3:5));
%             delta(delta>L-delta)=L-delta(delta>L-delta);
%             dist=[dist,(sum(delta.^2,2)).^0.5];
%             
% end 
% % H = zeros(64,49);
% % H1 = zeros(64,49,49);
% % 
% [s,~]=size(k);
% Q = sum((dist<cutoff),2)/s;
% % Q = reshape(Q,64,[]);
% nbins = 50;
% binedges = linspace(0,1,50);
% q_normalized1 = reshape(q_normalized',[],1);
% [e,f] = find(q_normalized1==0.15);
% figure
% h=histogram(Q(e),binedges,'Normalization','probability');
% plot(binedges(2:end),-log(h.Values),'--b')
% 
% % l = 1;
% % for i = 1:64
% %  h = histogram(Q(i,:)',binedges, 'Normalization', 'probability');
% %  H(i,:) = h.Values; 
% %  h1 = histogram2(q_normalized1,Q,'XBinEdges',binedges, 'YBinEDges', binedges, 'Normalization', 'probability');
% %  H1(i,:,:) = h1.Values; 
% %  l = l+64;
% % end
% % P = mean(H,1);
% % P1 = mean(H1,1);
% % P3(:,:) = P1(1,:,:);
% % plot(binedges(2:end),P)
% %  figure
% %   s=surf(binedges(2:end),binedges(2:end),-log(h1.Values)');
% %   s.EdgeColor = 'None';
% 
% 
% % dist = [];
% % 
% % parfor id1=1:length(k)
% %     id1
% %             delta = abs(trj_unfolded{k(id1,1)}(:,3:5)-trj_unfolded{k(id1,2)}(:,3:5));
% %             delta(delta>L-delta)=L-delta(delta>L-delta);
% %             dist=[dist,(sum(delta.^2,2)).^0.5];
% %             
% % end 
% % % H = zeros(max(n_unfolded),49);
% % [s,~]=size(k);
% % Qu = sum((dist<cutoff),2)/s;
% % ybinedges = linspace(0,1,25);
% % figure
% % l = 1;
% % nn = nan(max(n_unfolded),10001);
% % for i = 1:length(n_unfolded)
% %     nn(1:n_unfolded(i),i) = Qu(l:l+n_unfolded(i)-1);
% %      l = l+n_unfolded(i);
% % end 
% 
% % for i = 1:max(n_unfolded)
% %  h = histogram(nn(i,:)',ybinedges, 'Normalization', 'probability');
% %  H(i,:) = h.Values; 
% % %  l = l+n_unfolded(i);
% % % end
% % % P = mean(H,1);
% % plot(ybinedges(2:end),P)
% % h = histogram(Qu,ybinedges, 'Normalization', 'probability');
% % h.FaceColor = 'None';
% % h.EdgeColor = 'r';
% % h.LineStyle = '-.';
% % h.LineWidth = 2;
% 
% % dist = [];
% % 
% % parfor id1=1:length(k)
% %     id1
% %             delta = abs(trj_folded{k(id1,1)}(:,3:5)-trj_folded{k(id1,2)}(:,3:5));
% %             delta(delta>L-delta)=L-delta(delta>L-delta);
% %             dist=[dist,(sum(delta.^2,2)).^0.5];
% %             
% % end 
% % % H = zeros(64,49);
% % % [s,~]=size(k);
% % % Qf = sum((dist<cutoff),2)/s;
% % % xbinedges = linspace(0,0.42,50);
% % % ybinedges = linspace(0,1,50);
% % % figure
% % % l = 1;
% % % for i = 1:length(n_folded)
% % %  h = histogram(Qf(l:l+n_folded(i)-1)',ybinegdes, 'Normalization', 'probability');
% % %  H(i,:) = h.Values; 
% % %  l = l+n_folded(i);
% % % end
% % % P = mean(H,1);
% % % plot(ybinedges(2:end),P)
% % H = zeros(max(n_folded),49);
% % [s,~]=size(k);
% % Qu = sum((dist<cutoff),2)/s;
% % ybinedges = linspace(0,1,50);
% % figure
% % l = 1;
% % nn = nan(max(n_folded),10001);
% % for i = 1:length(n_folded)
% %     nn(1:n_folded(i),i) = Qu(l:l+n_folded(i)-1);
% %      l = l+n_folded(i);
% % end 
% % 
% % for i = 1:max(n_folded)
% %  h = histogram(nn(i,:)',ybinedges, 'Normalization', 'probability');
% %  H(i,:) = h.Values; 
% %  l = l+n_unfolded(i);
% % end
% % P = mean(H,1);
% % plot(ybinedges(2:end),P)
% toc
