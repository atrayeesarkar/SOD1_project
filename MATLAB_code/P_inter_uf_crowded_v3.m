clear all;
tic
L = 146.8;%in angstrom;
cutoff = 15;%in angstrom;
trj = cell(1,110);
contacts = 360;
% Loading the trjectories for 110 residues
parfor i = 1:110
    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_20/161/J_textfiles/','equil_',num2str(i),'.txt');
    trj{i}= load(t);
    trj{i}(640001:end,:)=[];
end 
s3 = 10000+1;
q= load('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_20/161/Q.dat');
q_n = q(s3:end-1,:);
q_n = q_n/contacts; 
clear q;
toc
disp("file loading done")
% size of trajectoriess
n=64;
q_normalized=q_n;
[t,~] = size(q_normalized); 
trj_unfolded = cell(1,110);
trj_folded = cell(1,110);
[i,j] = find(q_normalized<0.42);
k = [i,j];
k = sortrows(k,1);
n_unfolded = sum(q_normalized<0.42,2);
index1 = (k(:,1)-1)*64 +k(:,2);
[o,p] = find(q_normalized>0.42);
r = [o,p];
r = sortrows(r,1);
n_folded = sum(q_normalized>0.42,2);
index2 = (r(:,1)-1)*64 +r(:,2);
parfor l= 1:110
trj_unfolded{l} = trj{l}(index1,:);
trj_folded{l} = trj{l}(index2,:);
end
clear q_n q_normalized trj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%calcuating no of folded/unfolded proteins in each time frame , which is
%ssame for all cells in trj_unfolded or folded for the same time frame
%computing J(1,1)
[s_unfolded,~] = size(trj_unfolded{1}); 
[s_folded,~] = size(trj_folded{1}); 
Jaccard = zeros(110,t,110);
Jaccard1 = Jaccard;
J12_folded=[];
J12_unfolded=[];
J11_folded = [];
J11_unfolded = [];
%S = n_folded.*(n_folded-1);
%DD1 = zeros(109,sum(S));
%Dii = zeros(110,sum(S));
toc

parfor id1=1:109
	tic
   	 frame_no = 0;
    	 AA = zeros(1,t);
   	 i = 1;
%	 l = 1;
%	 Dij = zeros(109-id1+1,sum(S));
%    	 dii =zeros(1,sum(S));
    while i <s_folded
        frame_no = frame_no+1;
        index1 = nchoosek(i:i+n_folded(frame_no)-1,2);
      	index = [index1;index1(:,2),index1(:,1)];
        index = sortrows(index);
	delta = abs(trj_folded{id1}(index(:,1),3:5)-trj_folded{id1}(index(:,2),3:5)); 
	delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=(sum(delta.^2,2)).^0.5;
%        dii(l:l+S(frame_no)-1) = dist;
	intersectA=sum(dist<cutoff); % dist<cutoff gives 1 or 0  
        AA(frame_no)=intersectA;
        i = i+n_folded(frame_no);
%	l = l+S(frame_no);
    end
%    Dii(id1,:)=dii;
    J11_folded(id1,:)=AA;
     id1
     for id2 = id1+1:110
        frame_no = 0;
        AB = zeros(1,t);
        i = 1;
%	l = 1;
%	dij = zeros(1,sum(S));
        while i <s_folded
       	frame_no = frame_no+1;
	index1 = nchoosek(i:i+n_folded(frame_no)-1,2);
	index = [index1;index1(:,2),index1(:,1)];
       	index2 = sortrows(index);
        delta = abs(trj_folded{id1}(index2(:,1),3:5)-trj_folded{id2}(index2(:,2),3:5)); 
       	delta(delta>L-delta)=L-delta(delta>L-delta);
       	dist=(sum(delta.^2,2)).^0.5;
%       	dij(l:l+S(frame_no)-1) = dist';
	intersectAB=sum(dist<cutoff);
       	 AB(frame_no)=intersectAB;
    	i = i+n_folded(frame_no);
%       	l = l+S(frame_no);
	 end                                      
%	 Dij(id2,:) = dij;
       	 J12_folded =[J12_folded;AB];
    	end
%	DD=min(Dij,[],1);
%	DD1(id1,:)=DD;
toc   
end  
for id1=110:110
    frame_no = 0;
    AA = zeros(1,t);
    i = 1;
%    l = 1;
%    dii =zeros(1,sum(S));
    while i <s_folded
        frame_no = frame_no+1;
        index1 = nchoosek(i:i+n_folded(frame_no)-1,2);
        index = [index1;index1(:,2),index1(:,1)];
        index = sortrows(index);
	delta = abs(trj_folded{id1}(index(:,1),3:5)-trj_folded{id1}(index(:,2),3:5)); 
        delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=(sum(delta.^2,2)).^0.5;
%       dii(l:l+S(frame_no)-1) = dist';
	 intersectA=sum(dist<cutoff);
        AA(frame_no)=intersectA;
        i = i+n_folded(frame_no);
%        l = l+S(frame_no);
	end 
%     Dii(id1,:) = dii;
    J11_folded(id1,:)=AA;
end
%D1 = min(Dii,[],1);
%D2 = min(DD1,[],1);
%D3 = [D1;D2];
%D = min(D3,[],1);
%factor = zeros(1,t);
%l = 1;
%for i = 1:t
%    factor(i) = sum(D(l:l+S(i)-1)<cutoff,2);
%    l = l+S(i);
%end
k=0;
for id1=1:110
    Jaccard(id1,:,id1)=J11_folded(id1,:);
    if id1==110
        break;
    end
    for id2= id1+1:110
        k=k+1;
        Jaccard(id1,:,id2)=J12_folded(k,:);
    end 
end 
disp("all terms done");
for i=1:110
    for j = 1:110
    Jaccard(j,:,i)=Jaccard(i,:,j);
    end
end
disp('matrix1 done') 
disp('normalization done')
for i =1: length(n_folded)
	norm(i)=nchoosek(n_folded(i),2);
end


avg_J = (sum(Jaccard(:,:,:),2)/(2*sum(norm)));
disp('avg J done')
f1 = fopen('/uhpc/cheung/asarkar4/SOD_project/slurm_jobs/Jaccard_index/Pij_inter_crowded_folded-folded_WT_charged_161_15.txt_v6_modifiedtimeavg_withoutRmin','w');
fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
fclose(f1);
disp('all done')
toc
