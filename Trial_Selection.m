load('High_Data.mat','data_decision','data_outcome')

clearvars -except data_decision data_outcome pnums

st{1}=[5 5];
st{2}=[6 6];
% st{4}=[5 6];
% st{5}=[4 3];

pnums=[39 42 43 45 46 48 50 51 53 54 55 57 58 59 61 62 69 70 72 73 74 75 76 77 78 79 80];

load('electrodes_regions6.mat')

for i=1:27
els{i}=zeros(size(data_decision{i}.label,1),1);
for j=1:6
els{i}(elec_selec{j,pnums(i)})=j;
end
end

ps=27;
count=0;
for i=1:ps
    
count=count+1;

tr=[];
tr(:,1)=data_decision{i}.trialinfo(:,8)==1;
tr(:,2)=data_decision{i}.trialinfo(:,8)==2;
tr(:,3)=data_decision{i}.trialinfo(:,8)>2;
tr(:,4)=data_decision{i}.trialinfo(:,8)<3;

tr(:,5)=data_decision{i}.trialinfo(:,4)==1;
tr(:,6)=data_decision{i}.trialinfo(:,4)==2;
tr(:,7)=data_decision{i}.trialinfo(:,4)>=3;
tr(:,8)=data_decision{i}.trialinfo(:,4)<3;
tr(:,9)=data_decision{i}.trialinfo(:,4)>0;

% tr(:,7)=data_decision{i}.trialinfo(:,4)==1 & data_decision{i}.trialinfo(:,8)==1;
% tr(:,8)=data_decision{i}.trialinfo(:,4)==2 & data_decision{i}.trialinfo(:,8)==2;
% tr(:,9)=tr(:,7)==1 & tr(:,8)==1;

trialnumbers(i,:)=sum(tr,1);

for g=1:size(st,2)
for k=1:1

cfg=[];
cfg.channel=find(els{i}>0);
cfg.trials=find(tr(:,st{g}(k))==1);
mint=min(trialnumbers(i,[st{g}]));
cfg.trials=cfg.trials(1:mint);
cfg.keeptrials='yes';
trial_nums2(g,k,count)=mint;
decision_conditions{g,k,count}=ft_timelockanalysis(cfg,data_decision{i});
decision_outcome{g,k,count}=ft_timelockanalysis(cfg,data_outcome{i});

% for j=1:size(decision_outcome{g,k,count}.label,1)
% decision_outcome{g,k,count}.trial(:,j,:)=decision_outcome{g,k,count}.trial(:,j,:)-repmat(mean(decision_conditions{g,k,count}.trial(:,j,1:500),3),1,1,2501);
% decision_conditions{g,k,count}.trial(:,j,:)=decision_conditions{g,k,count}.trial(:,j,:)-repmat(mean(decision_conditions{g,k,count}.trial(:,j,1:500),3),1,1,2501);
% end

end
end
end