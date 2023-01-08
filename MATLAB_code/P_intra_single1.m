% % clear all;
tic %gives run time 
L = 140.0;%in angstrom;Length of 1 dimension of the simulation box
cutoff = 15;%in angstrom;% residues with distace less than cutoff are considered to be interacting
trj = cell(1,110); %initializing with cell data structure in matlab.It has 1 row of  110 cells. Each Cell can contain a matrix  
contacts=360;% total no. of native contacts

%%----------------------------------------WT_c------------------%
%
%% Loading the trjectories for 110 residues
%parfor i = 1:110 %parfor parallaelizes the for loop 
%    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_0/157_1/J_textfiles/','equil_',num2str(i),'.txt');%append is used because t has a differnt string value for different i value
%    trj{i}= load(t); %loading the trajectories inside the cell data structure t. Each cell in trj corresponds to one resddidue.There are 110 cells corresponding to 110 residues in total. Each cell contains a marix with rows corresponding to timeframes. The last 3 columns correspond to x,y,z coordinate of the protein for the residue (residue index is same as the cell index as mentioned earlier). The second column also gives the residue index.      
%      % trj{i}=trj{i}(1:320,:); %this was used as a test. Not part of the final code
%end
%t = []; 
%toc %gives run time so far
%q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_0/157_1/traj_comp_pbcmolcenter_WT_c.xtc.CA.Q");%loads q values (not normalized) 
%[s1,~]=size(trj{1}); % gives # timeframes. Note that trj contains only equilibiated time frames. 
%[s2,~]=size(q1);  % gives total no. of time frames including the non equilibriated ones.
%s3=s2-s1+1; % to discard the initial time frames before equilibrium
%q=q1(s3:s2); % to discard the initial time frames before equilibrium
%q_normalized = q/contacts; %normalizing q
%q1 = [];
%q = [];
%
%unfolded_time_frame = find(q_normalized<0.42); % indeces for unfolded time frames  
%folded_time_frame = find(q_normalized>0.42); %indices of folded time frames
%parfor i = 1:110 
%	trj_unfolded1{i}=trj{i}(unfolded_time_frame,:); % sorting out trajectories corresponding to unfolded and folded protein configurations.
%	trj_folded1{i}=trj{i}(folded_time_frame,:);
%end
%
%parfor i = 1:110 % same thing as above for run2 data
%    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_0/157_2/J_textfiles/','equil_',num2str(i),'.txt');
%    trj{i}= load(t);
%      % trj{i}=trj{i}(1:320,:);
%end
%t = [];
%toc
%q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_0/157_2/traj_comp_pbcmolcenter_WT_c.xtc.CA.Q");
%[s1,~]=size(trj{1});
%[s2,~]=size(q1);
%s3=s2-s1+1;
%q=q1(s3:s2);
%q_normalized = q/contacts;
%q1 = [];
%q = [];
%
%unfolded_time_frame = find(q_normalized<0.42);
%folded_time_frame = find(q_normalized>0.42);
%parfor i = 1:110
%        trj_unfolded2{i}=trj{i}(unfolded_time_frame,:);
%        trj_folded2{i}=trj{i}(folded_time_frame,:);
%end
%parfor i = 1:110 # same for run3 data
%    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_0/157_3/J_textfiles/','equil_',num2str(i),'.txt');
%    trj{i}= load(t);
%      % trj{i}=trj{i}(1:320,:);
%end
%t = [];
%toc
%q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/WT_charged/phic_0/157_3/traj_comp_pbcmolcenter_WT_c.xtc.CA.Q");
%[s1,~]=size(trj{1});
%[s2,~]=size(q1);
%s3=s2-s1+1;
%q=q1(s3:s2);
%q_normalized = q/contacts;
%q1 = [];
%q = [];
%unfolded_time_frame = find(q_normalized<0.42);
%folded_time_frame = find(q_normalized>0.42);
%parfor i = 1:110
%        trj_unfolded3{i}=trj{i}(unfolded_time_frame,:);
%        trj_folded3{i}=trj{i}(folded_time_frame,:);
%end
%
%parfor i = 1:110 % concatenated unfolded and folded trajectories from run1, run2 and run3 data
%	trj_unfolded{i}= [trj_unfolded1{i};trj_unfolded2{i};trj_unfolded3{i}];
%        trj_folded{i}=[trj_folded1{i};trj_folded2{i};trj_folded3{i}];
%end
%trj = [];
%trj_unfolded1 =[]; % freeing memory
%trj_unfolded2 =[];
%trj_unfolded3 =[];
%trj_folded1 =[];
%trj_folded2 =[];
%trj_folded3 =[];
%disp("file loading done")
%%---------------------WT_c folded------------------%
%% size of trajectoriess
%[s,~] = size(trj_folded{1}); %gives # timeframe
%%computing J(1,1)
%[t,~]=size(trj_folded{1}); % gives # timeframe 
%Jaccard=zeros(110,t,110); initializing
%J12=[];
%% k=0;
%parfor id1=1:109 % parallelizing for 
%     for id2 = id1+1:110
%        delta = abs(trj_folded{id1}(:,:3:5)-trj_folded{id2}(:,3:5)); % distance between residues 
%        delta(delta>L-delta)=L-delta(delta>L-delta); % periodic boundary condition
%        dist=(sum(delta.^2,2)).^0.5; 
%        AB=(dist(:)<cutoff);    
%        J12 =[J12;AB'];
%     end     
%end   
%    
%
%k=0;
%for id1=1:110
%    for id2= id1+1:110
%        k=k+1;             % error in this part of code
%        Jaccard(id1,:,id2)=J12(k,:);
%    end 
%end 
%disp("all terms done");
%
%for i=1:110
%    for j = 1:110
%    Jaccard(j,:,i)=Jaccard(i,:,j);
%    end
%end 
%delta = [];
%dist = [];
%AB = [];
%J12 = [];
%avg_J = sum(Jaccard(:,:,:),2)/t;
% f1 = fopen('Pij_WT_charged_157_15_intra_single_folded_v3.txt','w');
% fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
% fclose(f1);
%
%
%%--------------WT_c unfolded-------------------% 
%% size of trajectoriess
%[s,~] = size(trj_unfolded{1});
%%computing J(1,1)
%[t,~]=size(trj_unfolded{1});
%Jaccard=zeros(110,t,110);
%J12=[];
%% k=0;
%parfor id1=1:109
%     for id2 = id1+1:110
%        delta = abs(trj_unfolded{id1}(:,3:5)-trj_unfolded{id2}(:,3:5));
%        delta(delta>L-delta)=L-delta(delta>L-delta);
%        dist=(sum(delta.^2,2)).^0.5;
%        AB=(dist(:)<cutoff);
%        J12 =[J12;AB'];
%     end
%end
%
%
%k=0;
%for id1=1:110
%    for id2= id1+1:110
%        k=k+1;             % error in this part of code ?? fixed or not??
%        Jaccard(id1,:,id2)=J12(k,:);
%    end
%end
%disp("all terms done");
%
%for i=1:110
%    for j = 1:110
%    Jaccard(j,:,i)=Jaccard(i,:,j);
%    end
%end
%delta = [];
%dist = [];
%AB = [];
%J12 = [];
%
%avg_J = sum(Jaccard(:,:,:),2)/t;
% f1 = fopen('Pij_WT_charged_157_15_intra_single_unfolded_v3.txt','w');
% fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
% fclose(f1);
%

%-----------------------WT_u------------------%

% Loading the trjectories for 110 residues
parfor i = 1:110
    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_0/157_1/J_textfiles/','equil_',num2str(i),'.txt');
    trj{i}= load(t);
      % trj{i}=trj{i}(1:320,:);
end
t = [];
toc
q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_0/157_1/traj_comp_pbcmolcenter_WT_u.xtc.CA.Q");
[s1,~]=size(trj{1});
[s2,~]=size(q1);
s3=s2-s1+1;
q=q1(s3:s2);
q_normalized = q/contacts;
q1 = [];
q = [];

unfolded_time_frame = find(q_normalized<0.42);
folded_time_frame = find(q_normalized>0.42);
parfor i = 1:110
        trj_unfolded1{i}=trj{i}(unfolded_time_frame,:);
        trj_folded1{i}=trj{i}(folded_time_frame,:);
end

parfor i = 1:110
    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_0/157_2/J_textfiles/','equil_',num2str(i),'.txt');
    trj{i}= load(t);
      % trj{i}=trj{i}(1:320,:);
end
t = [];
toc
q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_0/157_2/traj_comp_pbcmolcenter_WT_u.xtc.CA.Q");
[s1,~]=size(trj{1});
[s2,~]=size(q1);
s3=s2-s1+1;
q=q1(s3:s2);
q_normalized = q/contacts;
q1 = [];
q = [];

unfolded_time_frame = find(q_normalized<0.42);
folded_time_frame = find(q_normalized>0.42);
parfor i = 1:110
        trj_unfolded2{i}=trj{i}(unfolded_time_frame,:);
        trj_folded2{i}=trj{i}(folded_time_frame,:);
end
parfor i = 1:110
    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_0/157_2/J_textfiles/','equil_',num2str(i),'.txt');
    trj{i}= load(t);
      % trj{i}=trj{i}(1:320,:);
end
t = [];
toc
q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/WT_uncharged/phic_0/157_2/traj_comp_pbcmolcenter_WT_u.xtc.CA.Q");
[s1,~]=size(trj{1});
[s2,~]=size(q1);
s3=s2-s1+1;
q=q1(s3:s2);
q_normalized = q/contacts;
q1 =[];
q = [];
unfolded_time_frame = find(q_normalized<0.42);
folded_time_frame = find(q_normalized>0.42);
parfor i = 1:110
        trj_unfolded3{i}=trj{i}(unfolded_time_frame,:);
        trj_folded3{i}=trj{i}(folded_time_frame,:);
end

parfor i = 1:110
        trj_unfolded{i}= [trj_unfolded1{i};trj_unfolded2{i};trj_unfolded3{i}];
        trj_folded{i}=[trj_folded1{i};trj_folded2{i};trj_folded3{i}];
end
trj = [];
trj_unfolded1 =[];
trj_unfolded2 =[];
trj_unfolded3 =[];
trj_folded1 =[];
trj_folded2 =[];
trj_folded3 =[];
disp("file loading done")
%%---------------------WT_u folded------------------%
%% size of trajectoriess
%[s,~] = size(trj_folded{1});
%%computing J(1,1)
%[t,~]=size(trj_folded{1});
%Jaccard=zeros(110,t,110);
%J12=[];
%% k=0;
%parfor id1=1:109
%     for id2 = id1+1:110
%        delta = abs(trj_folded{id1}(:,3:5)-trj_folded{id2}(:,3:5));
%        delta(delta>L-delta)=L-delta(delta>L-delta);
%        dist=(sum(delta.^2,2)).^0.5;
%        AB=(dist(:)<cutoff);
%        J12 =[J12;AB'];
%     end
%end
%
%
%k=0;
%for id1=1:110
%    for id2= id1+1:110
%        k=k+1;             % error in this part of code
%        Jaccard(id1,:,id2)=J12(k,:);
%    end
%end
%disp("all terms done");
%
%for i=1:110
%    for j = 1:110
%    Jaccard(j,:,i)=Jaccard(i,:,j);
%    end
%end
%delta = [];
%dist = [];
%AB = [];
%J12 = [];
%
%avg_J = sum(Jaccard(:,:,:),2)/t;
% f1 = fopen('Pij_WT_uncharged_157_15_intra_single_folded_v3.txt','w');
% fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
% fclose(f1);
%
%
%--------------WT_u unfolded-------------------% 
% size of trajectoriess
[s,~] = size(trj_unfolded{1});
%computing J(1,1)
[t,~]=size(trj_unfolded{1});
Jaccard=zeros(110,t,110);
J12=[];
% k=0;
parfor id1=1:109
     for id2 = id1+1:110
        delta = abs(trj_unfolded{id1}(:,3:5)-trj_unfolded{id2}(:,3:5));
        delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=(sum(delta.^2,2)).^0.5;
        AB=(dist(:)<cutoff);
        J12 =[J12;AB'];
     end
end


k=0;
for id1=1:110
    for id2= id1+1:110
        k=k+1;             % error in this part of code ?? fixed or not??
        Jaccard(id1,:,id2)=J12(k,:);
    end
end
disp("all terms done");

for i=1:110
    for j = 1:110
    Jaccard(j,:,i)=Jaccard(i,:,j);
    end
end

delta = [];
dist = [];
AB = [];
J12 = [];
avg_J = sum(Jaccard(:,:,:),2)/t;
 f1 = fopen('Pij_WT_uncharged_157_15_intra_single_unfolded_v3.txt','w');
 fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
 fclose(f1);


%----------------------------------------WT_m------------------%
%% Loading the trjectories for 110 residues
%parfor i = 1:110
%    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_1/J_textfiles/','equil_',num2str(i),'.txt');
%    trj{i}= load(t);
%      % trj{i}=trj{i}(1:320,:);
%end
%t = [];
%toc
%q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_1/traj_comp_pbcmolcenter_WT_m.xtc.CA.Q");
%[s1,~]=size(trj{1});
%[s2,~]=size(q1);
%s3=s2-s1+1;
%q=q1(s3:s2);
%parfor i = 1:110
%       trj{i}=trj{i}(1:9485,:);
%end
%q = q(2:9486);
%q_normalized = q/contacts;
%q1 = [];
%q = [];
%unfolded_time_frame = find(q_normalized<0.42);
%folded_time_frame = find(q_normalized>0.42);
%parfor i = 1:110
%        trj_unfolded{i}=trj{i}(unfolded_time_frame,:);
%        trj_folded{i}=trj{i}(folded_time_frame,:);
%end
%%
%%parfor i = 1:110
%%    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_2/J_textfiles/','equil_',num2str(i),'.txt');
%%    trj{i}= load(t);
%%      % trj{i}=trj{i}(1:320,:);
%%end
%%t = [];
%%toc
%%q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_2/traj_comp_pbcmolcenter_WT_m.xtc.CA.Q");
%%[s1,~]=size(trj{1});
%%[s2,~]=size(q1);
%%s3=s2-s1+1;
%%q=q1(s3:s2);
%%q_normalized = q/contacts;
%%q1 = [];
%%q = [];
%%unfolded_time_frame = find(q_normalized<0.42);
%%folded_time_frame = find(q_normalized>0.42);
%%parfor i = 1:110
%%        trj_unfolded2{i}=trj{i}(unfolded_time_frame,:);
%%        trj_folded2{i}=trj{i}(folded_time_frame,:);
%%end
%%
%%parfor i = 1:110
%%    t = append('/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_3/J_textfiles/','equil_',num2str(i),'.txt');
%%    trj{i}= load(t);
%%      % trj{i}=trj{i}(1:320,:);
%%end
%%t = [];
%%toc
%%q1 =load("/uhpc/cheung/asarkar4/SOD_project/metadata/G41D_mutant/phic_0/157_3/traj_comp_pbcmolcenter_WT_m.xtc.CA.Q");
%%[s1,~]=size(trj{1});
%%[s2,~]=size(q1);
%%s3=s2-s1+1;
%%q=q1(s3:s2);
%%q_normalized = q/contacts;
%%q1=[];
%%q = [];
%%unfolded_time_frame = find(q_normalized<0.42);
%%folded_time_frame = find(q_normalized>0.42);
%%parfor i = 1:110
%%        trj_unfolded3{i}=trj{i}(unfolded_time_frame,:);
%%        trj_folded3{i}=trj{i}(folded_time_frame,:);
%%end
%%
%%parfor i = 1:110
%%        trj_unfolded{i}= [trj_unfolded1{i};trj_unfolded2{i};trj_unfolded3{i}];
%%        trj_folded{i}=[trj_folded1{i};trj_folded2{i};trj_folded3{i}];
%%end
%trj = [];
%trj_unfolded1 =[];
%trj_unfolded2 =[];
%trj_unfolded3 =[];
%trj_folded1 =[];
%trj_folded2 =[];
%trj_folded3 =[];
%disp("file loading done")
%%---------------------WT_m folded------------------%
%% size of trajectoriess
%[s,~] = size(trj_folded{1});
%%computing J(1,1)
%[t,~]=size(trj_folded{1});
%Jaccard=zeros(110,t,110);
%J12=[];
%% k=0;
%parfor id1=1:109
%     for id2 = id1+1:110
%        delta = abs(trj_folded{id1}(:,3:5)-trj_folded{id2}(:,3:5));
%        delta(delta>L-delta)=L-delta(delta>L-delta);
%        dist=(sum(delta.^2,2)).^0.5;
%        AB=(dist(:)<cutoff);
%        J12 =[J12;AB'];
%     end
%end
%
%
%k=0;
%for id1=1:110
%    for id2= id1+1:110
%        k=k+1;             % error in this part of code
%        Jaccard(id1,:,id2)=J12(k,:);
%    end
%end
%disp("all terms done");
%
%for i=1:110
%    for j = 1:110
%    Jaccard(j,:,i)=Jaccard(i,:,j);
%    end
%end
%delta = [];
%dist = [];
%AB = [];
%J12 = [];
%avg_J = sum(Jaccard(:,:,:),2)/t;
% f1 = fopen('Pij_G41D_mutant_157_15_intra_single_folded_v3_157_1testQ_first9485frames.txt','w');
% fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
% fclose(f1);
%
%%%
%%--------------WT_m unfolded-------------------% 
%% size of trajectoriess
%[s,~] = size(trj_unfolded{1});
%%computing J(1,1)
%[t,~]=size(trj_unfolded{1});
%Jaccard=zeros(110,t,110);
%J12=[];
%% k=0;
%parfor id1=1:109
%     for id2 = id1+1:110
%        delta = abs(trj_unfolded{id1}(:,3:5)-trj_unfolded{id2}(:,3:5));
%        delta(delta>L-delta)=L-delta(delta>L-delta);
%        dist=(sum(delta.^2,2)).^0.5;
%        AB=(dist(:)<cutoff);
%        J12 =[J12;AB'];
%     end
%end
%
%
%k=0;
%for id1=1:110
%    for id2= id1+1:110
%        k=k+1;             % error in this part of code ?? fixed or not??
%        Jaccard(id1,:,id2)=J12(k,:);
%    end
%end
%disp("all terms done");
%
%for i=1:110
%    for j = 1:110
%    Jaccard(j,:,i)=Jaccard(i,:,j);
%    end
%end
%
%delta = [];
%dist = [];
%AB = [];
%J12 = [];
%avg_J = sum(Jaccard(:,:,:),2)/t;
% f1 = fopen('Pij_G41D_mutant_157_15_intra_single_unfolded_v3_157_1testQ_first9485frames.txt','w');
% fprintf(f1,'%1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f %1.10f\n',avg_J);
% fclose(f1);

%
%
%
%
%
%
%
%

toc

