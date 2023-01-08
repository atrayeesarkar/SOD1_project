%This is the latest version (05/30/2021) not the J1.m in the mac desktop
clear all;
tic
L = 146.8;%in angstrom;
cutoff = 15;%in angstrom;
trj1 = cell(1,110);
contacts = 360;
% Loading the trjectories for 110 residues
parfor i = 1:110
    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_20/161/J_textfiles/','equil_',num2str(i),'.txt');
    trj{i}= load(t);
    trj{i}(640001:end,:)=[];
%    trj{i}=trj{i}(1:64,:);
%trj{i}=trj{i}(1:512000,:);
end 

s3 = 10000 +1;

q= load('/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_20/161/Q.dat');
q_n = q(s3:end-1,:);
%q_n = q(s3,:);
q_n = q_n/contacts; 
%q_n = q_n(1:8000,:);
toc
disp("file loading done")
% size of trajectoriess
n=64;

q_normalized=q_n;
[t,~] = size(q_normalized); 

trj_unfolded = cell(1,110);
trj_folded = cell(1,110);
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
clear q_n q_normalized trj q
disp("traj sorting done")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%calcuating no of folded/unfolded proteins in each time frame , which is
%ssame for all cells in trj_unfolded or folded for the same time frame
%computing J(1,1)
[s_unfolded,~] = size(trj_unfolded{1}); 
[s_folded,~] = size(trj_folded{1}); 
% [~,t]=size(sample);
Jaccard_unfolded=zeros(110,t,110); % t = # of time frmes, assuming every timeframe has atleast one unfolded and unfolfded proteins
J12_unfolded=[];
Jaccard_folded=zeros(110,t,110); % t = # of time frmes, assuming every timeframe has atleast one unfolded and unfolfded proteins
J12_folded=[];
J11_unfolded = [];
J11_folded = [];
%Dii = [];
%Dij=[];
%S = n_folded.*(n_folded-1);
S = 2*n_unfolded.*n_folded;
%DD1 = zeros(109,sum(S));
%Dii = zeros(110,sum(S));
parfor id1=1:109
	tic
    frame_no = 0;
    AA = zeros(1,t);
    i = 1;
    j = 1;
%    l=1;
%    Dij = zeros(109-id1+1,sum(S));
    
%    dii = zeros(1,sum(S));
%    dii = [];
    while i <s_unfolded
        frame_no = frame_no+1;
%        a = [i,i:+n_unfolded(frame_no)-1];
%        b = [j,j:+n_folded(frame_no)-1];
	a = i:i+n_unfolded(frame_no)-1;
        b = j:j+n_folded(frame_no)-1;
        [x,y] = ndgrid(a,b);
        index= [x(:),y(:)];
        index1 = sortrows(index);
        index2 = sortrows(index,2);
        delta1 = abs(trj_unfolded{id1}(index1(:,1),3:5)-trj_folded{id1}(index1(:,2),3:5));
        delta2 = abs(trj_folded{id1}(index2(:,2),3:5)-trj_unfolded{id1}(index2(:,1),3:5));
        delta = [delta1;delta2];
%        index(any(diff(sort(index,2),[],2)==0,2),:)=[];
%        [~,idx] = unique(sort(index,2),'rows','stable');
%        index = index(idx,:);
%         index = sortrows(index);
%        index = nchoosek(i:i+n_unfolded(frame_no)-1,2);
%        delta = abs(trj_unfolded{id1}(index(:,1),3:5)-trj_folded{id1}(index(:,2),3:5)); 
        delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=(sum(delta.^2,2)).^0.5;
%        dii = [dii,dist'];
%	dii(l:l+S(frame_no)-1) = dist';
        intersectA=sum(dist<cutoff); % dist<cutoff gives 1 or 0  
        AA(frame_no)=intersectA;
        i = i+n_unfolded(frame_no);
%	l = l+S(frame_no);
        j = j+n_folded(frame_no);
    end
    J11_unfolded(id1,:)=AA;
%    Dii(id1,:)= dii;	
%clear a b x y index1 index2 delta1 deta2 delta dist dii intersectA AA 
%dii = [];
intersectA = [];
AA = [];
     id1
     for id2 = id1+1:110
        frame_no = 0;
        AB = zeros(1,t);
        i = 1;
	j = 1;
        l = 1;
%        dij = zeros(1,sum(S));
%        dij = [];
        while i <s_unfolded
        frame_no = frame_no+1;
%         index2 = unique(nchoosek([i:i+n_unfolded(frame_no)-1,i:i+n_unfolded(frame_no)-1],2),"rows");
%         index2(any(diff(sort(index2,2),[],2)==0,2),:)=[];
%        a = [i,i:+n_unfolded(frame_no)-1];
%        b = [j,j:+n_folded(frame_no)-1];
	a = i:i+n_unfolded(frame_no)-1;
        b = j:j+n_folded(frame_no)-1;
        [x,y] = ndgrid(a,b);
        index= [x(:),y(:)];
	index1 = sortrows(index);
        index2 = sortrows(index,2);
        delta1 = abs(trj_unfolded{id1}(index1(:,1),3:5)-trj_folded{id2}(index1(:,2),3:5));
        delta2 = abs(trj_folded{id1}(index2(:,2),3:5)-trj_unfolded{id2}(index2(:,1),3:5));
        delta = [delta1;delta2];
%        index2(any(diff(sort(index2,2),[],2)==0,2),:)=[];
%	index2 = sortrows(index2);
%        delta = abs(trj_unfolded{id1}(index2(:,1),3:5)-trj_folded{id2}(index2(:,2),3:5)); 
        delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=(sum(delta.^2,2)).^0.5;
%	dij = [dij,dist'];
%	dij(l:l+S(frame_no)-1) = dist';
        intersectAB=sum(dist<cutoff);
        AB(frame_no)=intersectAB;
        i = i+n_unfolded(frame_no);
	 j = j+n_folded(frame_no);
%	 l = l+S(frame_no);
        end  
        J12_unfolded =[J12_unfolded;AB];
%	Dij = [Dij,dij];
%	Dij(id2,:) = dij;
%clear a b x y index1 index2 delta1 deta2 delta dist dij intersectAB AB
%dij = [];
intersectAB = [];
AB = [];
	  end   
%	DD=min(Dij,[],1);
%      DD1(id1,:)=DD;
%clear Dij DD
%Dij = [];
%DD = [];
toc
delete(gcp('nocreate'))
end  
    
for id1=110:110
    frame_no = 0;
    AA = zeros(1,t);
    i = 1;
    j = 1;
%    l = 1;
%    dii = zeros(1,sum(S));
%    dii = [];
    while i <s_unfolded
        frame_no = frame_no+1;
%        a = [i,i:+n_unfolded(frame_no)-1];
%        b = [j,j:+n_folded(frame_no)-1];
	a = i:i+n_unfolded(frame_no)-1;
        b = j:j+n_folded(frame_no)-1;
        [x,y] = ndgrid(a,b);
        index= [x(:),y(:)];
%        [~,idx] = unique(sort(index,2),'rows','stable');
%        index = index(idx,:);
%	index = sortrows(index);
%        index = nchoosek(i:i+n_unfolded(frame_no)-1,2);
%        delta = abs(trj_unfolded{id1}(index(:,1),3:5)-trj_folded{id1}(index(:,2),3:5)); 
	index1 = sortrows(index);
        index2 = sortrows(index,2);
        delta1 = abs(trj_unfolded{id1}(index1(:,1),3:5)-trj_folded{id1}(index1(:,2),3:5));
        delta2 = abs(trj_folded{id1}(index2(:,2),3:5)-trj_unfolded{id1}(index2(:,1),3:5));
        delta = [delta1;delta2];
        delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=(sum(delta.^2,2)).^0.5;
%        dii = [dii,dist'];
%	dii(l:l+S(frame_no)-1) = dist';
        intersectA=sum(dist<cutoff); % dist<cutoff gives 1 or 0  
        AA(frame_no)=intersectA;
        i = i+n_unfolded(frame_no);
        j = j+n_folded(frame_no);
%	l = l+S(frame_no);
    end
%    Dii(id1,:)=dii;
    J11_unfolded(id1,:)=AA;
end
clear a b x y index1 index2 delta1 deta2 delta dist  intersectA AA

%[S,~] = size(dist);
%D1 = min(Dii,[],1);
%D2 = min(DD1,[],1);
%D3 = [D1;D2];
%D = min(D3,[],1);
%factor = sum(reshape(D<cutoff,[S,t]),1);
%factor = zeros(1,t);
%l = 1;
%for i = 1:t
%    factor(i) = sum(D(l:l+S(i)-1)<cutoff,2);
%    l = l+S(i);
%end
k=0;
for id1=1:110
    Jaccard(id1,:,id1)=J11_unfolded(id1,:);
    if id1==110
        break;
    end
    for id2= id1+1:110
        k=k+1;
        Jaccard(id1,:,id2)=J12_unfolded(k,:);
    end 
end 
disp("all terms done");

for i=1:110
    for j = 1:110
    Jaccard(j,:,i)=Jaccard(i,:,j);
    end
end 
disp("mm started")
%mm=repmat(max(Jaccard,[],[1,3]),110,1,110); %total no of proteins in contact in each frame
%mm = max(Jaccard,[],[1,3]);
disp('mm done')
%for i = 1:t
%        Jaccard1(:,i,:) = Jaccard(:,i,:)/mm(i);
%	Jaccard1(:,i,:) = Jaccard(:,i,:)/factor(i);
%end
disp('normalizatio done done')
%Jaccard1=vpa(Jaccard./mm);
for i =1: length(n_unfolded)
	norm(i)= n_unfolded(i)*n_folded(i);
end
1
avg_J = (sum(Jaccard(:,:,:),2)/(2*sum(norm)));


2
f1 = fopen('/uhpc/cheung/asarkar4/SOD_project/slurm_jobs/Jaccard_index/Pij_inter_crowded_unfolded-folded_G41D_mutant_161_15.txt_v5_withoutRmin','w');
fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
fclose(f1);


%f2 = fopen('Dii_WT_u_inter_uf.txt','w');
%for i = 1:size(Dii,1)
%    fprintf(f2,'%5.8f\t',Dii(i,:));
%    fprintf(f2,'\n');
%end
%fclose(f2)
%f3 = fopen('Dij_WT_u_inter_uf.txt','w');
%for i = 1:size(Dij,1)
%    fprintf(f3,'%5.8f\t',Dij(i,:));
%    fprintf(f3,'\n');
%end
%fclose(f3)
%
%
%
%
%f2=fopen('unfolded_WT_u_inter_uf.txt','w');
%for i = 1:size(kk,1)
%    fprintf(f2,'%5.8f\t',kk(i,:));
%    fprintf(f2,'\n');
%end
%fclose(f2)
%f3 = fopen('folded_WT_u_inter_uf.txt','w');
%for i = 1:size(r,1)
%    fprintf(f3,'%5.8f\t',r(i,:));
%    fprintf(f3,'\n');
%end
%fclose(f3)
toc

