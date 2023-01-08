% % % clear all;
% tic
% L = 146.8;%in angstrom;
% cutoff = 15;%in angstrom;
% trj = cell(1,110);
% contacts = 360;
% % Loading the trjectories for 110 residues
% parfor i = 1:110
%     t = append('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_20/161_1/J_textfiles/','equil_',num2str(i),'.txt');
%     trj{i}= load(t);
%     trj{i}(640001:end,:)=[];
% end
% s3 = 10000+1;
% q= load('/Volumes/TAB_RESEARCH_CHEUNG/Atrayee/SOD_project/Data/G41D_Mutant/phic_20/161_1/Q.dat');
% q_n = q(s3:end-1,:);
% q_n = q_n/contacts;
% clear q;
% toc
% disp("file loading done")
% % size of trajectoriess
% n=64;
% q_normalized=q_n;
% [t,~] = size(q_normalized);
% trj_unfolded = cell(1,110);
% trj_folded = cell(1,110);
% [i,j] = find(q_normalized<0.42);
% k = [i,j];
% k = sortrows(k,1);
% n_unfolded = sum(q_normalized<0.42,2);
% index1 = (k(:,1)-1)*64 +k(:,2);
% [o,p] = find(q_normalized>0.42);
% r = [o,p];
% r = sortrows(r,1);
% n_folded = sum(q_normalized>0.42,2);
% index2 = (r(:,1)-1)*64 +r(:,2);
% parfor l= 1:110
% trj_unfolded{l} = trj{l}(index1,:);
% trj_folded{l} = trj{l}(index2,:);
% end
% Indu = index1;
% Indf = index2;

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
id1=41;
     for id2 =1:8
        frame_no = 0;
        AB = zeros(1,t);
        i = 1;
        j = 1;
        while i <s_unfolded
        frame_no = frame_no+1;
        a = i:i+n_unfolded(frame_no)-1;
        b = j:j+n_folded(frame_no)-1;
        [x,y] = ndgrid(a,b);
        index= [x(:),y(:)];
        index1 = sortrows(index);
        delta = abs(trj_unfolded{id1}(index1(:,1),3:5)-trj_folded{id2}(index1(:,2),3:5));
        delta(delta>L-delta)=L-delta(delta>L-delta);
        dist=(sum(delta.^2,2)).^0.5;
        I = find(dist<cutoff);
        F1{frame_no}=index1(I,:);
        F1{frame_no}(:,1) = k(F1{frame_no}(:,1),2); %r in place of k for folded
        F1{frame_no}(:,2) = r(F1{frame_no}(:,2),2);
        i = i+n_unfolded(frame_no);
         j = j+n_folded(frame_no);
        end
        
        R{id2}= F1;
     end 
toc






        