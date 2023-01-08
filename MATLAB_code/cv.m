id = 170;
aa = [140,150,155,156,157,158,159,160,170,180];
e = cell(1,length(aa));
sp_ht_u = zeros(1,length(aa));
sp_ht_m = zeros(1,length(aa));
sp_ht_c = zeros(1,length(aa));
for i = 1:length(id)
    [l,k]= find(aa == id(i));
    t = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/','',num2str(id(i)),'_1/','energy.txt');
    e{k}= load(t);
    e{k} = e{k}(10001:end,:);
end 


for i = 1:length(id)
    [l,k]= find(aa == id(i));
    j = var(e{k}(:,2));
    sp_ht_u(k) = j;
end

for i = 1:length(id)
    [l,k]= find(aa == id(i));
    t = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/','',num2str(id(i)),'_1/','energy.txt');
    e{k}= load(t);
    e{k} = e{k}(10001:end,:);
end 

for i = 1:length(id)
    [l,k]= find(aa == id(i));
    j = var(e{k}(:,2));
    sp_ht_c(k) = j;
    
end

 for i = 1:length(id)
    [l,k]= find(aa == id(i));
    t = append('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/','',num2str(id(i)),'_1/','energy.txt');
    e{k}= load(t);
    e{k} = e{k}(10001:end,:);
end 


for i = 1:length(id)
    [l,k]= find(aa == id(i));
    j = var(e{k}(:,2));
    sp_ht_m(k) = j
    
end



id = [140,150,155,156,157,158,159,160,180];
for i = 1:length(id)
    [l,k]= find(aa == id(i));
    t1 = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/','',num2str(id(i)),'_1/','energy.txt');
    t2 = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/','',num2str(id(i)),'_2/','energy.txt');
    t3 = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/','',num2str(id(i)),'_3/','energy.txt');
    e1 = load(t1);
    e1 = e1(10001:end,:);
    e2 = load(t2);
    e2 = e2(10001:end,:);
    e3 = load(t3);
    e3 = e3(10001:end,:);
    e{k}= [e1;e2;e3];
    
end 

for i = 1:length(id)
    [l,k]= find(aa == id(i));
    j = var(e{k}(:,2));
    sp_ht_u(k) = j;
end
figure
subplot(2,1,1); 
axis square;
hold on
grid on
plot(aa,sp_ht_u)
for i = 1:length(id)
    [l,k]= find(aa == id(i));
    t1 = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/','',num2str(id(i)),'_1/','energy.txt');
    t2 = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/','',num2str(id(i)),'_2/','energy.txt');
    t3 = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/','',num2str(id(i)),'_3/','energy.txt');
    e1 = load(t1);
    e1 = e1(10001:end,:);
    e2 = load(t2);
    e2 = e2(10001:end,:);
    e3 = load(t3);
    e3 = e3(10001:end,:);
    e{k}= [e1;e2;e3];
    
end 


for i = 1:length(id)
    [l,k]= find(aa == id(i));
    j = var(e{k}(:,2));
    sp_ht_c(k) = j;
    
end


 for i = 1:length(id)
    [l,k]= find(aa == id(i));
    t1 = append('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/','',num2str(id(i)),'_1/','energy.txt');
    t2 = append('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/','',num2str(id(i)),'_2/','energy.txt');
    t3 = append('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/','',num2str(id(i)),'_3/','energy.txt');
    e1 = load(t1);
    e1 = e1(10001:end,:);
    e2 = load(t2);
    e2 = e2(10001:end,:);
    e3 = load(t3);
    e3 = e3(10001:end,:);
    e{k}= [e1;e2;e3];
    
end 

for i = 1:length(id)
    [l,k]= find(aa == id(i));
    j = var(e{k}(:,2));
    sp_ht_m(k) = j
    
end
plot(aa,sp_ht_m) 
 plot(aa,sp_ht_c) 
%  %%%%%%%%%%% crowded cases %%%%%%%%%%%%%
%  
%  id = [120,130,140,145,150,155,158,159,160,161,162,163,164,165,170,180,190,200,210,220];
% % id1 = [-40,-20,-10,-5,-4,-3,-2,-1,0,10,20,30,40,50,60,20,80,90,100, id, 300, 310, 320,330,340,350,360,370,380,390,400,401,402,403,404,405, 410, 420, 440];
% 
% parfor i = 1:20
%     t = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_20/','',num2str(id(i)),'_1/','energy.txt');
%     e{i}= load(t);
%     e{i} = e{i}(10001:end,:);
% end 
% sp_ht_u = [];
% k = [];
% for i = 1:20
%     j = var(e{i}(:,2));
%     sp_ht_u = [sp_ht_u,j];
% end
%  k = repelem(id,round(sp_ht_u));
% %  figure
% subplot(2,1,2); 
% hold on
% grid on
% [f,x1,b] = ksdensity(k,'bandwidth',20);
% %  plot(x1,f,'--b')
% plot(id,sp_ht_u)
%  axis 'square';
% parfor i = 1:20
%     t = append('/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_20/','',num2str(id(i)),'_1/','energy.txt');
%     e{i}= load(t);
% end 
% sp_ht_c = [];
% k = [];
% for i = 1:20
%     j = var(e{i}(:,2));
%     sp_ht_c = [sp_ht_c,j];
%     
% end
% k = repelem(id,round(sp_ht_c));
% [f,x1,b] = ksdensity(k,'bandwidth',20);
% %  plot(x1,f,': k')
% plot(id,sp_ht_c) 
%  parfor i = 1:20
%     t = append('/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_20/','',num2str(id(i)),'_1/','energy.txt');
%     e{i}= load(t);
%     
% end 
% sp_ht_m = [];
% k = [];
% for i = 1:20
%     j = var(e{i}(:,2));
%     sp_ht_m = [sp_ht_m,j];
%     
% end
% k = repelem(id,round(sp_ht_c));
% [f,x1,b] = ksdensity(k,'bandwidth',20);
% %  plot(x1,f,'-. k')
% plot(id,sp_ht_m)
%  axis 'square';
%  legend ('WT_u','WT_c','WT_m')