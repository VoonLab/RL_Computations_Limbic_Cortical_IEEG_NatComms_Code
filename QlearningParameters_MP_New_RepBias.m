% Parameter estimation for Q-learning model
% Church & Harrison modification of Pessiglione

% data structure
% There is 1 .mat file per session, hence 3 per subjects, containing a cell
% array named 'data' (clever naming) with the following columns :
% data(:,1) = session [1:3]
% data(:,2) = trial [1:90]
% data(:,3) = pair [1=gain; 2=look; 3=loss]
% data(:,4) = onset time
% data(:,5) = onset slice
% data(:,6) = response [1=go; -1=nogo]
% data(:,7) = choice [1=correct; -1=incorrect]
% data(:,8) = feedback [1=good; -1=bad]
% data(:,9) = reaction time [0 is for nogo]
% data(:,10) = hesitation [release after press]

clearvars

pnums=[35 39 42 43 44 45 46 48 49 50 51 53 54 55 56 57 58 59 61 62 69 70 72 73 74 75 76 77 78 79 80 81 82 84 86];

prob2=zeros(150,86);

proba1=zeros(150,86);
probb1=zeros(150,86);
qu3=zeros(150,86);
qa3=zeros(150,86);
qb3=zeros(150,86);

for ii=1:length(pnums)

pnums=[35 39 42 43 44 45 46 48 49 50 51 53 54 55 56 57 58 59 61 62 69 70 72 73 74 75 76 77 78 79 80 81 82 84 86];
    
%LLB=zeros(10,10); % initialise vectors for the loglikelihood of kappa x lambda
final = [];
behav_filename=strcat('/Volumes/Luis_HDD_3/SEEG_All/Patients/Patient_',num2str(pnums(ii)),'/behav_data/PLT',num2str(pnums(ii)),'.mat');
load(behav_filename);
outcomes=outcomes';
gain_loss=gain_loss';
side=side';
resp_cat=resp_cat';
Acc1=zeros(150,1);
Acc1(gain_loss==1 & outcomes==1)=1;
%Acc1(gain_loss==2 & outcomes~=2)=1;
%Acc1(gain_loss==1 & outcomes~=1)=-1;
Acc1(gain_loss==2 & outcomes==2)=-1;
Acc2=zeros(150,1);
Acc2(gain_loss==1 & resp_cat==side)=1;
Acc2(gain_loss==1 & resp_cat~=side)=-1;
Acc2(gain_loss==2 & resp_cat==side)=-1;
Acc2(gain_loss==2 & resp_cat~=side)=1;
sess=[ones(75,1);ones(75,1)+1];
data=[gain_loss Acc2 Acc1 sess zeros(150,1)+pnums(ii) outcomes];

osel=[1 -1];

PE=zeros(150,1);
q2=zeros(150,1);
qu2=zeros(150,1);
qa2=zeros(150,1);
qb2=zeros(150,1);

for jj=1:2 %reward loss
    for kk=1:2 %sess
        
        selectdata=data(data(:,1)==jj & data(:,4)==kk,:);
        inx=find(data(:,1)==jj & data(:,4)==kk);
        LL=zeros(100,100,100);  % initialise vectors for the loglikelihood of alpha x beta

% loop through all combinations of alpha and beta
for i=1:100
    alpha=i/100;  % thus 0 <= alpha <= 1
    for j=1:100
        beta=j/100; % thus 0 <= beta <= 1
        for k=1:100
            theta=k/100;
            disp([ii jj kk i j k])

        qA=0;  % intialise q values to zero
        qB=0;

        proba=[];  % record the probabilities corresponding to the chosen actions
        error=[]; % record prediction errors (outcome - expectation)
        
        pchose=0;

        % loop through trials
        for t=1:length(selectdata)

            % calculate probabilities of choice A and choice B using softmax function
        if pchose==1
        pA=exp(qA+theta/beta)/(exp(qA+theta/beta)+exp(qB/beta));
        pB=exp(qB/beta)/(exp(qA+theta/beta)+exp(qB/beta));
        end

        if pchose==2
        pA=exp(qA/beta)/(exp(qA/beta)+exp(qB+theta/beta));
        pB=exp(qB+theta/beta)/(exp(qA/beta)+exp(qB+theta/beta));
        end

        if pchose==0
        pA=exp(qA/beta)/(exp(qA/beta)+exp(qB/beta));
        pB=exp(qB/beta)/(exp(qA/beta)+exp(qB/beta));
        end

            if jj==1
            reward=selectdata(t,3)==1;
            end

            if jj==2
            reward=-(selectdata(t,3)==-1);
            end

            if selectdata(t,2)==1 % correct choice
                proba(t)=log(pA);  % note p(chosen action)
                % update q values, arbitrary reward value set at 1
                error(t)=reward-qA;
                qA=qA+alpha*error(t);
                pchose=1;

            else % incorrect choice, chose B
                proba(t)=log(pB); % note p(chosen action)
                % update q values, arbitrary reward value set at 1
                error(t)=reward-qB;
                qB=qB+alpha*error(t);
                pchose=2;
                    
            end

        end

        % update likelihhod array
        LL(i,j,k)=LL(i,j,k)+sum(proba);
        end

    end

end

% [alpha,beta]=find(LL==max(max(LL))); % find optimal values of alpha and beta
% alpha=alpha/100; % thus 0 <= alpha <= 1
% beta=beta/100; % thus 0 <= beta <= 1

p1=LL(1,1,1);
coords1=[1 1 1];
for i=1:100
    for j=1:100
        for k=1:100
        p2=LL(i,j,k);
        if p2>p1
            p1=LL(i,j,k);
            coords1=[i j k];
        end
        end
    end
end

alpha=coords1(1)/100;
beta=coords1(2)/100;
theta=coords1(3)/100;

LLF(pnums(ii),jj,kk)=LL(coords1(1),coords1(2),coords1(3));

alpha2(pnums(ii),jj,kk)=alpha;
beta2(pnums(ii),jj,kk)=beta;
theta2(pnums(ii),jj,kk)=theta;

clearvars -except theta2 theta LLF proba1 probb1 osel jj kk ii alpha beta selectdata PE2 rew_loss data inx PE data2 PE3 pnums q2 q Q3 trials PE4 Q4 prob1 prob2 p2 alpha2 beta2 qu qu3 qu2 qa2 qa3 qb2 qb3 qa3_4 qb3_4 qu3_4

qA=0;  % intialise q values to zero
qB=0;

proba=[];  % record the probabilities corresponding to the chosen actions
error=[]; % record prediction errors (outcome - expectation)
p2=[];
pchose=0;

% loop through trials
for t=1:length(selectdata)

    % calculate probabilities of choice A and choice B using softmax function
        if pchose==1
        pA=exp(qA+theta/beta)/(exp(qA+theta/beta)+exp(qB/beta));
        pB=exp(qB/beta)/(exp(qA+theta/beta)+exp(qB/beta));
        end

        if pchose==2
        pA=exp(qA/beta)/(exp(qA/beta)+exp(qB+theta/beta));
        pB=exp(qB+theta/beta)/(exp(qA/beta)+exp(qB+theta/beta));
        end

        if pchose==0
        pA=exp(qA/beta)/(exp(qA/beta)+exp(qB/beta));
        pB=exp(qB/beta)/(exp(qA/beta)+exp(qB/beta));
        end
    
    if jj==1
    reward=selectdata(t,3)==1;
    end
    
    if jj==2
    reward=-(selectdata(t,3)==-1);
    end
    
    proba12(t)=pA;
    probb12(t)=pB;

    if selectdata(t,2)==1 % correct choice, chose A
        proba(t)=log(pA);  % note p(chosen action)
        p2(t)=pA;
        % update q values, arbitrary reward value set at 1
        error(t)=reward-qA;
        qA=qA+alpha*error(t);
        q(t)=qA;
        qu(t)=qB;
        pchose=1;

    else % incorrect choice, chose B
        proba(t)=log(pB); % note p(chosen action)
        p2(t)=pB;
        % update q values, arbitrary reward value set at 1
        error(t)=reward-qB;
        qB=qB+alpha*error(t);
        q(t)=qB;
        qu(t)=qA;
        pchose=2;

    end
    
    qa1(t)=qA;
    qb1(t)=qB;

end

PE(inx)=error;
q2(inx)=q;
prob2(inx,pnums(ii))=p2;
proba1(inx,pnums(ii))=proba12;
probb1(inx,pnums(ii))=probb12;
qu2(inx)=qu;

qa2(inx)=qa1;
qb2(inx)=qb1;

end
end

PE2(:,ii)=[PE];
PE3(:,pnums(ii))=PE;
Q3(:,pnums(ii))=q2;
qu3(:,pnums(ii))=qu2;
qa3(:,pnums(ii))=qa2;
qb3(:,pnums(ii))=qb2;

data2(:,:,pnums(ii))=data;

PE4(1:150,pnums(ii))=0;
Q4(1:150,pnums(ii))=0;

qa3_4(1:150,pnums(ii))=0;
qb3_4(1:150,pnums(ii))=0;
qu3_4(1:150,pnums(ii))=0;

for i=1:6
    clear PET PE_n1 g QC Q_n1 qa3_1 qa3_2 qb3_1 qb3_2
        g=find(trials==i);
        
        PE_n1=PE3(:,pnums(ii));
        PET=PE_n1(g);
        
        Q_n1=Q3(:,pnums(ii));
        QC=Q_n1(g);
        
        qa3_1=qa3(:,pnums(ii));
        qa3_2=qa3_1(g);
        
        qb3_1=qb3(:,pnums(ii));
        qb3_2=qb3_1(g);
        
        qu3_1=qu3(:,pnums(ii));
        qu3_2=qu3_1(g);
        
        PE4(g(2:end),pnums(ii))=[PET(1:end-1)];
        Q4(g(2:end),pnums(ii))=[QC(1:end-1)];
        
        qa3_4(g(2:end),pnums(ii))=[qa3_2(1:end-1)];
        qb3_4(g(2:end),pnums(ii))=[qb3_2(1:end-1)];
        qu3_4(g(2:end),pnums(ii))=[qu3_2(1:end-1)];
end

clearvars -except theta2 LLF jj kk ii PE2 data2 PE3 Q3 PE4 trials Q4 prob2 prob1 osel alpha2 beta2 proba1 probb1 qu3 qa3 qb3 qa3_4 qb3_4 qu3_4

end

save('PredictionErrors_rp.mat','PE3','Q3','PE4','Q4','prob2','alpha2','beta2','proba1','probb1','qu3','qa3','qb3','qa3_4','qb3_4','qu3_4','data2','LLF','theta2')

% for i=1:16
%     figure
%     subplot(1,2,1)
%     plot(PE2(squeeze(data2(:,1,i))==1,i));
%     title(num2str(i))
%     subplot(1,2,2)
%     plot(PE2(squeeze(data2(:,1,i))==2,i));
%     pause
%     close(gcf)
% end
