function [Distance, FundID_Least, FundID_Most] = FcmGrnn_NonFoF_2SampleTtestA(netattack,FundID,N)

% This function not only compute the similiarity distance among funds but
% also identify least & most distressed funds which will be used to conduct two-sample difference t-tests.
% The last column of netattack is the clasifier of event: 1=Alive,2=Liquidate
P1=netattack;
T1=P1(:,end)';
P1(:,end)=[];

[R1,~]=size(P1);
csum=20;  

data=P1;
[~,U,~] = fcm(data,2);    
for i=1:R1
    [~,idx]=max(U(:,i));
    a1(i)=idx;
end

Confusion_Matrix_FCM=zeros(3,3);
Confusion_Matrix_FCM(1,:)=[0:2];
Confusion_Matrix_FCM(:,1)=[0:2]';
for nf=1:2
    for nc=1:2
        Confusion_Matrix_FCM(nf+1,nc+1)=length(find(a1(find(T1==nf))==nc));
    end
end

cent1=P1(find(a1==1),:);cent1=mean(cent1);
cent2=P1(find(a1==2),:);cent2=mean(cent2);

for n=1:R1
    ecent1(n)=norm(P1(n,:)-cent1);
    ecent2(n)=norm(P1(n,:)-cent2);
end

for n=1:csum
    [va me1]=min(ecent1);
    [va me2]=min(ecent2);
    ecnt1(n,:)=P1(me1(1),:);ecent1(me1(1))=[];tcl(n)=1;
    ecnt2(n,:)=P1(me2(1),:);ecent2(me2(1))=[];tc2(n)=2;
end
P2=[ecnt1;ecnt2];T2=[tcl,tc2];
k=0;

LoopNum = N;
for nit=1:LoopNum
    net = newgrnn(P2',T2,4.5);   
    a2=sim(net,P1') ; 
    A2=a2;
    a2(find(a2<1.5))=1;
    a2(find(a2>=1.5))=2;
    cent1=P1(find(a2==1),:);cent1=mean(cent1);
    cent2=P1(find(a2==2),:);cent2=mean(cent2);
    
    ID_L = FundID(find(a2==1)); % Least distressed funds
    ID_M = FundID(find(a2==2)); % Most distressed funds
    
    for n=1:R1
        ecent1(n)=norm(P1(n,:)-cent1);
        ecent2(n)=norm(P1(n,:)-cent2);
    end
    
    if nit == LoopNum
        Distance(1:2,:) = [ecent1; ecent2];
    end
      
    for n=1:csum
        [va me1]=min(ecent1);
        [va me2]=min(ecent2);
        ecnt1(n,:)=P1(me1(1),:);ecent1(me1(1))=[];tc1(n)=1;
        ecnt2(n,:)=P1(me2(1),:);ecent2(me2(1))=[];tc2(n)=2;
        % Identifying those leaset and most similar fundID
        FundID_Least(n,:) = FundID(me1(1));
        FundID_Most(n,:)  = FundID(me2(1));
    end
    
    P2=[ecnt1;ecnt2];T2=[tc1,tc2];

    Confusion_Matrix_GRNN=zeros(3,3);
    Confusion_Matrix_GRNN(1,:)=[0:2];
    Confusion_Matrix_GRNN(:,1)=[0:2]';
    for nf=1:3
        for nc=1:3
            Confusion_Matrix_GRNN(nf+1,nc+1)=length(find(a2(find(T1==nf))==nc));
        end
    end
end


