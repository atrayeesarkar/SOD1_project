PUu = load('Pij_inter_crowded_unfolded-unfolded_WT_uncharged_161_15.txt_v6_modifiedtimeavg_withoutRmin');
PCu = load('Pij_inter_crowded_unfolded-unfolded_WT_charged_161_15.txt_v6_modifiedtimeavg_withoutRmin');

c1 = zeros(110,110);
[i,j]=find(PUu>1.25*PCu);
k = [i,j];
k1 = k;
for i = 1:length(k)
    c1(k(i,1),k(i,2))=40; % deep red
end 
[i,j]=find(PUu>1.10*PCu & PUu<1.25*PCu) ;
k = [i,j];
for i = 1:length(k)
    c1(k(i,1),k(i,2))=35;
end 

[i,j]=find(PCu>1.5*PUu);

k = [i,j];
for i = 1:length(k)
    c1(k(i,1),k(i,2))=5; % deep blue
end 
[i,j]=find(PCu>1.10*PUu & PCu<1.10*PUu);
k = [i,j];
for i = 1:length(k)
    c1(k(i,1),k(i,2))=15;
end 


PUu = load('Pij_inter_crowded_unfolded-folded_WT_uncharged_161_15.txt_v5_withoutRmin');
PCu = load('Pij_inter_crowded_unfolded-folded_WT_charged_161_15.txt_v5_withoutRmin');

c2 = zeros(110,110);
[i,j]=find(PUu>1.25*PCu);
k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=40; % deep red
end 
[i,j]=find(PUu>1.10*PCu & PUu<1.25*PCu) ;
k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=35;
end 

[i,j]=find(PCu>1.25*PUu);

k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=5; % deep blue
end 
[i,j]=find(PCu>1.10*PUu & PCu<1.25*PUu);
k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=15;
end 

PUu = load('Pij_inter_crowded_folded-folded_WT_uncharged_161_15.txt_v6_modifiedtimeavg_withoutRmin');
PCu = load('Pij_inter_crowded_folded-folded_WT_charged_161_15.txt_v6_modifiedtimeavg_withoutRmin');

c3 = zeros(110,110);
[i,j]=find(PUu>1.25*PCu);
k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=40; % deep red
end 
[i,j]=find(PUu>1.10*PCu & PUu<1.25*PCu) ;
k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=35;
end 

[i,j]=find(PCu>1.25*PUu);

k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=5; % deep blue
end 
[i,j]=find(PCu>1.10*PUu & PCu<1.25*PUu);
k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=15;
end 



figure
h=heatmap(triu(reshape(c1,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
h.NodeChildren(3).YDir='normal';  
h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
% cmap=[winter(500);white(1);flipud(autumn(499))];
% colormap(c1);
figure
h=heatmap(triu(reshape(c2,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
h.NodeChildren(3).YDir='normal';  
h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
% cmap=[winter(500);white(1);flipud(autumn(499))];
% colormap(c2);
figure
h=heatmap(triu(reshape(c3,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
h.NodeChildren(3).YDir='normal';  
h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
% cmap=[winter(500);white(1);flipud(autumn(499))];
% colormap(c3);






PGu = load('Pij_inter_crowded_unfolded-unfolded_G41D_mutant_161_15.txt_v6_modifiedtimeavg_withoutRmin');
PCu = load('Pij_inter_crowded_unfolded-unfolded_WT_charged_161_15.txt_v6_modifiedtimeavg_withoutRmin');

c1 = zeros(110,110);
[i,j]=find(PGu>1.25*PCu);
k = [i,j];
Kr = k;
for i = 1:length(k)
    c1(k(i,1),k(i,2))=40; % deep red
end 
[i,j]=find(PGu>1.10*PCu & PGu<1.25*PCu) ;
k = [i,j];
for i = 1:length(k)
    c1(k(i,1),k(i,2))=35;
end 

[i,j]=find(PCu>1.25*PGu);

k = [i,j];
Kb = k;
for i = 1:length(k)
    c1(k(i,1),k(i,2))=5; % deep blue
end 
[i,j]=find(PCu>1.10*PGu & PCu<1.25*PGu);
k = [i,j];
for i = 1:length(k)
    c1(k(i,1),k(i,2))=15;
end 


PGu = load('Pij_inter_crowded_unfolded-folded_G41D_mutant_161_15.txt_v5_withoutRmin');
PCu = load('Pij_inter_crowded_unfolded-folded_WT_charged_161_15.txt_v5_withoutRmin');

c2 = zeros(110,110);
[i,j]=find(PGu>1.25*PCu);
k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=40; % deep red
end 
[i,j]=find(PGu>1.10*PCu & PGu<1.25*PCu) ;
k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=35;
end 

[i,j]=find(PCu>1.25*PGu);

k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=5; % deep blue
end 
[i,j]=find(PCu>1.10*PGu & PCu<1.25*PGu);
k = [i,j];
for i = 1:length(k)
    c2(k(i,1),k(i,2))=15;
end 

PGu = load('Pij_inter_crowded_folded-folded_G41D_mutant_161_15.txt_v6_modifiedtimeavg_withoutRmin');
PCu = load('Pij_inter_crowded_folded-folded_WT_charged_161_15.txt_v6_modifiedtimeavg_withoutRmin');

c3 = zeros(110,110);
[i,j]=find(PGu>1.25*PCu);
k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=40; % deep red
end 
[i,j]=find(PGu>1.10*PCu & PGu<1.25*PCu) ;
k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=35;
end 

[i,j]=find(PCu>1.25*PGu);

k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=5; % deep blue 
end 
[i,j]=find(PCu>1.10*PGu & PCu<1.25*PGu);
k = [i,j];
for i = 1:length(k)
    c3(k(i,1),k(i,2))=15;
end 



figure
h=heatmap(triu(reshape(c1,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
h.NodeChildren(3).YDir='normal';  
h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
% cmap=[winter(500);white(1);flipud(autumn(499))];
% colormap(c1);
figure
h=heatmap(triu(reshape(c2,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
h.NodeChildren(3).YDir='normal';  
h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
% cmap=[winter(500);white(1);flipud(autumn(499))];
% colormap(c2);
figure
h=heatmap(triu(reshape(c3,110,110)),'GridVisible','off','MissingDataColor',[1 1 1]);
h.NodeChildren(3).YDir='normal';  
h.XDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
h.YDisplayLabels = {'','','','','','','','','','','','','','','','','','','','20','','','','','','','','','','','','','','','','','','','','40','','','','','','','','','','','','','','','','','','','','60','','','','','','','','','','','','','','','','','','','', '80','','','','','','','','','','','','','','','','','','','', '100','','','','','','','','','',''};
% cmap=[winter(500);white(1);flipud(autumn(499))];
% colormap(c3);


