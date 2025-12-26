%========================================================================================================================%
%    This code is to calcualte the distance (Fuzzy Cluster Means) between Fund-of-Funds and other trading style funds every year
%     and combine other explanatory variables ready to be for SAS saving data for aML estimation                      
%========================================================================================================================%
CHAR = xlsread('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\HedgeFundCharacters.xlsx','FundCharacteristics');
     %-- Creating dummies for Style, LegalStructure	ManagementFeePayablePeriod	SubscriptionFrequency	RedemptionFrequency	DomicileState
StyleD = dummyvar(CHAR(:,2));
LegalD = dummyvar(CHAR(:,27));
MnpfeD = dummyvar(CHAR(:,28));    
SubscD = dummyvar(CHAR(:,29));
RedemD = dummyvar(CHAR(:,30));
StateD = dummyvar(CHAR(:,31));
         %-- Need to group trading style into five categories
StyleDR = StyleD(:,2)+StyleD(:,3)+StyleD(:,8)+StyleD(:,9);    % Directional
StyleRV = StyleD(:,1)+StyleD(:,4)+StyleD(:,6);                % Relative Value
StyleED = StyleD(:,5)+StyleD(:,10)+StyleD(:,11)+StyleD(:,12); % Event_Driven
StyleFF = StyleD(:,7);                                        % Fund of Fund
StyleNN = StyleD(:,13)+StyleD(:,14);                          % Others + Not Defind

LegalNN = LegalD(:,45);                                        % Not defined
LegalOE = LegalD(:,46)+LegalD(:,47)+LegalD(:,48);              % OEID
LegalLI = LegalD(:,39)+LegalD(:,40)+LegalD(:,41);              % Limited Parterner/Corporation/Liability
LegalEM = LegalD(:,19)+LegalD(:,20)+LegalD(:,21)+LegalD(:,22); % Exempted
LegalD(:,[19:22 39:41 45:48])=[];
LegalOT = sum(LegalD,2);                                       % Others

MnpfeNN = MnpfeD(:,7);                                          % Not defined
MnpfeMM = MnpfeD(:,6);                                          % Monthly
MnpfeQQ = MnpfeD(:,8);                                          % Quarterly
MnpfeYY = MnpfeD(:,1);                                          % Annually
MnpfeD(:, [1 6 7 8])=[];
MnpfeOT = sum(MnpfeD,2);                                        % Others

SubscMM = SubscD(:,5);                                          % Monthly
SubscNN = SubscD(:,6);                                          % Not defined
SubscD(:, [5 6])=[];
SubscOT = sum(SubscD,2);                                        % Others

RedemMM = RedemD(:,7);                                          % Monthly
RedemFN = RedemD(:,6);                                          % Fornightly
RedemNN = RedemD(:,8);                                          % Not defined
RedemD(:,[6 7 8])=[];
RedemOT = sum(RedemD,2);                                        % Others

StateNN = StateD(:,31);                                         % Not defined
StateDL = StateD(:,8);                                          % Delaware
StateD(:,[8 31])=[];
StateOT = sum(StateD,2);                                        % Other states

CHAR(:,[2:4 16 27:32])=[];
CHAR = [CHAR StyleDR StyleRV StyleED StyleFF StyleNN LegalNN LegalOE LegalLI LegalEM LegalOT MnpfeNN ...
             MnpfeMM MnpfeQQ MnpfeYY MnpfeOT SubscMM SubscNN SubscOT RedemMM RedemFN RedemNN RedemOT StateNN StateDL StateOT]; 
FundStyle = [CHAR(:,1) StyleFF];                                %[FundID Fund_of_Fund];
%=== Hedge Fund Time Series Performance ===%
DATA_TS = xlsread('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\HedgeFundPerformanceFlow.xlsx');
% DATA_TS = [FundID	Year Month Day Event Return	NAV	AUM	FundFlow Flow_YM Flow_YM_SP	Flow_YM_DP]
     %-- Cleaning data: Deleting funds less than 24 months --%
b=1; counter=0;
while b<=size(DATA_TS,1)
    counter = counter + 1;
    n = length(find(DATA_TS(:,1) == DATA_TS(b,1)));
    LifeSpan(b:b+n-1,1) = n;
    b = b + n;
end
DATA_TS(LifeSpan<=24,:) = [];                % Funds less than 24 months
LifeSpan(LifeSpan<=24,:) = [];
Statistics.LifeSpan = [mean(LifeSpan) std(LifeSpan) median(LifeSpan) min(LifeSpan) max(LifeSpan)]; % Structure
% AFTER CLEANING, THERE ARE NOW ??? FUNDS ALTOGETHER INCLUDING  ??? LIQUIDATED AND ??? ALIVE FUNDS%
%-- Splitting funds into Non_FOF funds and FOF funds --%
for i = 1:size(DATA_TS,1)
    fof(i,1) = FundStyle(ismember(FundStyle(:,1),DATA_TS(i,1)),2);
end
DATA_IN = DATA_TS(fof==0,:);      % Non-FoF funds
DATA_OUT = DATA_TS(fof==1,:);     % FoF Funds

%== Selecting both alive and liquidated hedge funds for every year starting from 1993 to 2021 ===%
T{1} = DATA_IN(DATA_IN(:,2)<1993,:);       % Those Non-FOF funds before 1993
Tfof{1} = DATA_OUT(DATA_OUT(:,2)<1993,:);  % Those FOF funds 
Year = 1993:2021;
%W=vertcat([], T{1});
for t = 1:length(Year)
    T{t+1} = DATA_IN(DATA_IN(:,2)==Year(t),:);
   % W = vertcat(W,T{t+1});
    Tfof{t+1} = DATA_OUT(DATA_OUT(:,2)==Year(t),:);
end
%=== Compute Distance between every fund and the centriod for each cluster ===%
for k = 1:size(T,2)
    CT = CHAR(ismember(CHAR(:,1),unique(T{k}(:,1))),:);
    %--- Substituting missing characters by the average for each variable ---%
    for i = 1:size(CT,2)-1
        CT(CT(:,i+1)==-999,i+1) = mean(CT(CT(:,i+1)~=-999,i+1));
    end
    %--- Computing fund Return, Std, Size, Flow averages ---%
    b=1; counter = 0;
    while b<=size(T{k},1)
        counter = counter+1;
        n = length(find(T{k}(:,1) == T{k}(b,1)));
        if n>1
            performance(counter,:) = [mean(T{k}(b:b+n-1,6)) std((T{k}(b:b+n-1,6))) mean(T{k}(b:b+n-1,8:9))];
        else
            performance(counter,:) = [T{k}(b:b+n-1,6) 0 T{k}(b:b+n-1,8:9)];
        end
        classfier(counter,1) = T{k}(b+n-1,5) + 1; % 1 = Alive; 2 = Liquidated
        ID(counter,1) = T{k}(b,1);
        b = b + n;
    end
    X = [performance classfier];
    %  GUI_YI X
    X(:,1) =(X(:,1)-min(X(:,1)))/(max(X(:,1))-min(X(:,1)));
    X(:,2) =(X(:,2)-min(X(:,2)))/(max(X(:,2))-min(X(:,2)));
    X(:,3) =(X(:,3)-min(X(:,3)))/(max(X(:,3))-min(X(:,3)));
    X(:,4) =(X(:,4)-min(X(:,4)))/(max(X(:,4))-min(X(:,4)));
        
    % Preparing for FOF Dataset
    b=1; counter = 0;
    while b<=size(Tfof{k},1)
        counter = counter+1;
        n = length(find(Tfof{k}(:,1) == Tfof{k}(b,1)));
        if n>1
            Y(counter,:) = [mean(Tfof{k}(b:b+n-1,6)) std((Tfof{k}(b:b+n-1,6))) mean(Tfof{k}(b:b+n-1,8:9))];
        else
            Y(counter,:) = [Tfof{k}(b:b+n-1,6) 0 Tfof{k}(b:b+n-1,8:9)];
        end
        %class(counter,1) = Tfof{k}(b+n-1,5) + 1; % 1 = Alive; 2 = Liquidated
        IDfof(counter,1) = Tfof{k}(b,1);
        Pefof(counter,:) = Y(counter,:);
        b = b + n;
    end
    %  GUI_YI Y
    Y(:,1) =(Y(:,1)-min(Y(:,1)))/(max(Y(:,1))-min(Y(:,1)));
    Y(:,2) =(Y(:,2)-min(Y(:,2)))/(max(Y(:,2))-min(Y(:,2)));
    Y(:,3) =(Y(:,3)-min(Y(:,3)))/(max(Y(:,3))-min(Y(:,3)));
    Y(:,4) =(Y(:,4)-min(Y(:,4)))/(max(Y(:,4))-min(Y(:,4)));
    % Calculating the distance using FCM-GRNN algorithm
    Distance = FcmGrnn_HedgeFund_OutSample(X,Y);
    Distance_Result{k} = [IDfof Distance' Pefof]; % Distance for each year from 1994 to 2018 for each FOF
    CT=[]; performance=[]; classfier=[]; ID=[]; IDfof=[]; X=[]; Y=[]; Distance=[]; Pefof=[];
end
%Stack the datasets, container into one dataset
Result=[];
for k=1:length(Distance_Result)
    d=Distance_Result{k};
    d = [d(:,1) (1991+k)*ones(length(d),1) d(:,2:end)]; %[FundID Year Distance1 Distance2 Return STD AUM Flow];
    Result=padconcatenation(Result,d,1);
    d=[]; 
end
RD = sortrows(Result,[1,2]); % FCM-NeuralNet Algorithm Distance results between FOF with other trading strategy funds

%% Now organise financial market conditions(i.e., VIX, CreditSpread, IlliquidIndex) Time-Varying variables
TV = xlsread('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\MarketConditionVariables2022.xlsx');
% TV = [Year Month VIX BAA10YM	AggLiq InnovLiq TradedLiq CRSP]
Year=[];
Year = 1992:2021;
for i = 1:length(Year)
    for j = 1:6 %[VIX BAA10YM AggLiq InnovLiq TradedLiq CRSP]
        TimeVarying_A(i,j) = mean(TV(TV(:,1)==Year(i),j+2));
    end
end
TimeVarying = [Year' TimeVarying_A];
%--- Merge RD with TimeVarying variables
for i = 1:size(RD,1)
    RDTV(i,:) = [RD(i,:) TimeVarying(TimeVarying(:,1)==RD(i,2),2:end)];% Matching by Year
end
% RD_TV = [FundID Year Distance1 Distance2 Return STD AUM Flow VIX BAA10YM AggLiq InnovLiq TradedLiq CRSP] %

xlswrite('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\FoFData4SAS.xlsx', RDTV);


    

    








    
    
    
    
    











    
