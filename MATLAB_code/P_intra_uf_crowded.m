clear all;
tic
L = 146.8;%in angstrom;
cutoff = 15;%in angstrom;
trj = cell(1,110);
contacts = 360; % total no. of native contacts
% Loading the trjectories for 110 residues
parfor i = 1:110
  t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_20/161/J_textfiles/','equil_',num2str(i),'.txt');
  trj{i}= load(t);
 % trj{i}(end-63:end,:)=[];
  %trj{i}(end-63:end,:)=[];
   trj{i}=trj{i}(1:64,:);
end 
q=[];
q_n = [];
s3 = 10000 +1;
%for i = 1:64
%  Q = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_20/161/Q/','q_',num2str(i),'.dat');
% q(:,i))= load(Q);
q = load('/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_20/161/Q.dat');
%   q(1:s3,i)=[];
%  q_n(:,i) = q(s3:end,i);
%  q_n(:,i) = q_n(:,i)/contacts; 
%end
%q_n = q(s3:end-1,:);
q_n = q(s3,:);
  q_n = q_n./contacts; 
toc
disp("file loading done")
% size of trajectoriess
n=64;
%sample = [];
[s,~] = size(trj{1}) 
%u1 = trj{1}(1,1);
%u2 = trj{1}(s,1);
%for u = u1:874:u2
%sample = [sample,u];
%end
%q_normalized=q_n(sample-10000,:);
q_normalized = q_n;
[t,~] = size(q_normalized) 
%parfor i = 1:110
%    found = ismember(trj{i}(:,1),sample);
%    trj{i}= trj{i}(found,:)
%end
trj_unfolded = cell(1,110);
trj_folded = cell(1,110);
%n_unfolded=zeros(1,t);
%n_folded=zeros(1,t);
%for i = 1:t
%  k = (i-1)*n;
%  for j = 1:n
%    if q_normalized(i,j)<.42
%      n_unfolded(i)=n_unfolded(i)+1;
%      parfor l = 1:110
%        trj_unfolded{l} = [trj_unfolded{l};trj{l}(j+k,:)];
%      end 
%    elseif q_normalized(i,j)>.42
%      n_folded(i) = n_folded(i)+1;
%      parfor l = 1:110
%        trj_folded{l} = [trj_folded{l};trj{l}(j+k,:)];
%      end 
%    end 
%  end
%end
[i,j] = find(q_normalized<0.42);
k = [i,j];
k = sortrows(k,1);%sort matrix axccording to ascending oreder of column 1
n_unfolded = sum(q_normalized<0.42,2);
index1 = (k(:,1)-1)*64 +k(:,2);
[o,p] = find(q_normalized>0.42);
r = [o,p];
r = sortrows(r,1) ;
n_folded = sum(q_normalized>0.42,2);
index2 = (r(:,1)-1)*64 +r(:,2);
parfor l= 1:110
trj_unfolded{l} = trj{l}(index1,:);
trj_folded{l} = trj{l}(index2,:);
end
size(trj_unfolded{1})
size(trj_folded{1})
disp("traj sorting done") 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%calcuating no of folded/unfolded proteins in each time frame , which is
%ssame for all cells in trj_unfolded or folded for the same time frame
%computing J(1,1)
[s_unfolded,~] = size(trj_unfolded{1}); 
[s_folded,~] = size(trj_folded{1}); 
%[~,t]=size(sample);
Jaccard_unfolded=zeros(110,t,110); % t = # of time frmes, assuming every timeframe has atleast one folded and unfolfded proteins
J12_unfolded=[];
Jaccard_folded=zeros(110,t,110); % t = # of time frmes, assuming every timeframe has atleast one folded and unfolfded proteins
J12_folded=[];
% k=0;
parfor id1=1:109
    id1
   for id2 = id1+1:110
    delta_unfolded = abs(trj_unfolded{id1}(:,3:5)-trj_unfolded{id2}(:,3:5)); 
    delta_unfolded(delta_unfolded>L-delta_unfolded)=L-delta_unfolded(delta_unfolded>L-delta_unfolded);
    dist_unfolded=(sum(delta_unfolded.^2,2)).^0.5;
    delta_folded = abs(trj_folded{id1}(:,3:5)-trj_folded{id2}(:,3:5)); 
    delta_folded(delta_folded>L-delta_folded)=L-delta_folded(delta_folded>L-delta_folded);
    dist_folded=(sum(delta_folded.^2,2)).^0.5;
  i_unfolded=1;
  i_folded=1;
  frame_no=0;
  AB_unfolded=zeros(1,t);
  AB_folded=zeros(1,t);
  while i_unfolded< s_unfolded
    frame_no=frame_no+1;
    %AB_unfolded(frame_no)=sum(dist_unfolded(i_unfolded:i_unfolded+n_unfolded(frame_no)-1)<cutoff)/n_unfolded(frame_no);
    AB_unfolded(frame_no)=sum(dist_unfolded(i_unfolded:i_unfolded+n_unfolded(frame_no)-1)<cutoff);
    i_unfolded = i_unfolded+n_unfolded(frame_no);
  end   
  i=1;
  frame_no=0;
  while i_folded < s_folded
    frame_no=frame_no+1;
    %AB_folded(frame_no)=sum(dist_folded(i_folded:i_folded+n_folded(frame_no)-1)<cutoff)/n_folded(frame_no);
    AB_folded(frame_no)=sum(dist_folded(i_folded:i_folded+n_folded(frame_no)-1)<cutoff);
    i_folded = i_folded+n_folded(frame_no);
  end  
    J12_unfolded =[J12_unfolded;AB_unfolded];
    J12_folded =[J12_folded;AB_folded];
   end
end   
k=0;
for id1=1:110    %parfor loop won't work here
  for id2= id1+1:110
    k=k+1;
    Jaccard_unfolded(id1,:,id2)=J12_unfolded(k,:);
    Jaccard_folded(id1,:,id2)=J12_folded(k,:);
  end 
end 
disp("all  terms done");
for i=1:110      % parfor loop won't work here
  for j = 1:110
  Jaccard_unfolded(j,:,i)=Jaccard_unfolded(i,:,j);
  Jaccard_folded(j,:,i)=Jaccard_folded(i,:,j);
  end
end 
%Jaccard1_unfolded=vpa(Jaccard_unfolded./(64));
%Jaccard1_folded=vpa(Jaccard_folded./(64));

avg_J_unfolded = sum(Jaccard_unfolded(:,:,:),2)/(sum(n_unfolded));
avg_J_folded = sum(Jaccard_folded(:,:,:),2)/(sum(n_folded));
 f1 = fopen('Pij_intra_crowded_unfolded_G41D_mutant_161_15_new5_first1frame.txt','w');
 fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J_unfolded);
 fclose(f1);
 f1 = fopen('Pij_intra_crowded_folded_G41D_mutant_161_15_new5_first1frame.txt','w');
 fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J_folded);
 fclose(f1);
%figure

%h=heatmap(triu(reshape(avg_J_unfolded,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
%h.NodeChildren(3).YDir='normal';  
%h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
%h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
%figure
%h=heatmap(triu(reshape(avg_J_folded,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
%h.NodeChildren(3).YDir='normal';  
%h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
%h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
%savefig("WT_uncharged_160v2_15_sampling.fig")
%saveas(gcf,"WT_uncharged_161_15v2_sampling.fig")
toc
