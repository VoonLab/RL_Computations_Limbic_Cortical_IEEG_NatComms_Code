function [trl] = trialfun_plt(cfg)

load('plt_marks.mat')
load(behav_filename)
load('PredictionErrors_rp.mat')

r1=find(gain_loss==1);
r2=zeros(150,1);
r2(r1)=r1;
r3=zeros(150,1);
r3(r1(2:end))=r1(1:end-1);

l1=find(gain_loss==2);
l2=zeros(150,1);
l2(l1)=l1;
l3=zeros(150,1);
l3(l1(2:end))=l1(1:end-1);

u2=r2+l2;
u3=r3+l3;

Q3=Q3(:,partic);
PE3=PE3(:,partic);
PE4=PE4(:,partic);
Q4=Q4(:,partic);
prob2=prob2(:,partic);
proba1=proba1(:,partic);
probb1=probb1(:,partic);
qu3=qu3(:,partic);

qa3=qa3(:,partic);
qb3=qb3(:,partic);
qa3_4=qa3_4(:,partic);
qb3_4=qb3_4(:,partic);
qu3_4=qu3_4(:,partic);

Q3_2=Q3;
PE3_2=PE3;
PE4_2=PE4;
Q4_2=Q4;
prob2_2=prob2;
qu3_2=qu3;
qa3_4_2=qa3_4;
qb3_4_2=qb3_4;

% qa3_2=qa3;
% qb3_2=qb3;
% qa3_4_2=qa3_4;
% qb3_4_2=qb3_4;

% Q3_2(gain_loss==2 & outcomes==2,:)=0-Q3(gain_loss==2 & outcomes==2);
% PE3_2(gain_loss==2 & outcomes==2,:)=0-PE3(gain_loss==2 & outcomes==2);
% PE4_2(gain_loss==2 & outcomes==2,:)=0-PE4(gain_loss==2 & outcomes==2);
% Q4_2(gain_loss==2 & outcomes==2,:)=0-Q4(gain_loss==2 & outcomes==2);
% prob2_2(gain_loss==2 & outcomes==2,:)=1-prob2(gain_loss==2 & outcomes==2);

%Q3_2(gain_loss==2,:)=0-Q3(gain_loss==2); %important
%PE3_2(gain_loss==2,:)=0-PE3(gain_loss==2);

PE4_2(gain_loss==2,:)=0-PE4(gain_loss==2);
Q4_2(gain_loss==2,:)=0-Q4(gain_loss==2);
prob2_2(gain_loss==2,:)=1-prob2(gain_loss==2);
qu3_2(gain_loss==2,:)=0-qu3(gain_loss==2);
% qa3_4_2(gain_loss==2,:)=0-qa3_4(gain_loss==2);
% qb3_4_2(gain_loss==2,:)=0-qb3_4(gain_loss==2);

if partic>35
for i=1:6
rpt(find(trials==i),1)=1:25;
end
end

if partic==35
for i=1:3
rpt(find(trials==i),1)=1:50;
end
end

if partic~=58

cfg.dataset=filename2;
hdr=ft_read_header(cfg.dataset);
event=ft_read_event(cfg.dataset);

for j=1:length(event)
    sample(j)=event(j).sample;
    vals{j}=event(j).value;
end

vals(1)=[];
sample(1)=[];

sample=round(sample);

for j=1:length(vals)
    valuesv(j)=sscanf(vals{j},'S%d');
end

if partic==43
    valuesv(1)=[];
end

if partic==82
    valuesv(valuesv==61)=189;
    valuesv(valuesv==62)=190;
end

if partic==86
    valuesv(valuesv==185)=189;
    valuesv(valuesv==186)=190;
    valuesv(valuesv==27)=31;
    valuesv(valuesv==185)=967;
    valuesv(valuesv==186)=967;
end

f1=find(valuesv==189);
sample(f1-1)=[];
valuesv(f1-1)=[];
f1=find(valuesv==189);
sample(f1)=[];
valuesv(f1)=[];

f1=find(valuesv==190);
sample(f1-2)=[];
valuesv(f1-2)=[];
f1=find(valuesv==190);
sample(f1-1)=[];
valuesv(f1-1)=[];
f1=find(valuesv==190);
sample(f1)=[];
valuesv(f1)=[];

trlbegin(:,1)=sample(valuesv<=3)-4000;
trlend(:,1)=sample(valuesv>=31 & valuesv<=35)+6000;
offset=zeros(size(trlbegin(:,1),1),1)-4000;
rtsamp=sample(valuesv==21 | valuesv==22);
trl_nums2(1:150,1)=1:150;
patient(1:150,1)=partic;

if partic==86
rtsamp(:,1:150)=0;
end

elseif partic==58

load('Triggers.mat')

trigs=(double(int64(trigs))+1265019+10000)/10000;
[valuesv sample]=findpeaks(trigs);

f1=find(valuesv(2:end)==189)+1;
sample(f1-1)=[];
valuesv(f1-1)=[];
f1=find(valuesv(2:end)==189)+1;
sample(f1)=[];
valuesv(f1)=[];

valuesv(1)=[];
sample(1)=[];
sample=round((sample/2500)*1000);

f1=find(valuesv==190);
sample(f1-2)=[];
valuesv(f1-2)=[];
f1=find(valuesv==190);
sample(f1-1)=[];
valuesv(f1-1)=[];
f1=find(valuesv==190);
sample(f1)=[];
valuesv(f1)=[];

trlbegin(:,1)=sample(valuesv<=3)-4000;
trlend(:,1)=sample(valuesv>=31)+6000;
offset=zeros(size(trlbegin(:,1),1),1)-4000;
rtsamp=sample(valuesv==21 | valuesv==22);
trl_nums2(1:150,1)=1:150;
patient(1:150,1)=partic;

end

trials2=[0 trials(1:end-1)];

ouc=zeros(150,1);
ouc(gain_loss==1 & outcomes==1)=1;
ouc(gain_loss==2 & outcomes==3)=1;

ouc(gain_loss==1 & outcomes==3)=2;
ouc(gain_loss==2 & outcomes==2)=2;

sess=[ones(1,75) ones(1,75)+1];

ouc2=zeros(150,1);
for ik=1:3
    for ij=1:2
        tt=find(gain_loss==ik & sess==ij);
        ouc2(tt(2:end))=ouc(tt(1:end-1));
    end
end

trl=[trlbegin trlend offset trlbegin trlend rtsamp' gain_loss' side' rt' resp_cat' outcomes' ouc2 trl_nums2 trials' trials2' total_money(:,1:3) patient qu3_2 qa3 qb3 qa3_4_2 qb3_4_2 qu3_4 proba1 probb1 prob2 prob2_2 u2 u3 PE3 Q3 PE4 Q4 PE3_2 Q3_2 PE4_2 Q4_2 rpt];

clearvars -except trl event