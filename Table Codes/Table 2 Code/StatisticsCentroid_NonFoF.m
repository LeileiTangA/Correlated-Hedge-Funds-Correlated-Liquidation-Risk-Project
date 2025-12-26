%==========================================================================================================%
%    This code is to test statistic difference between distressed and successful Non-FoFfunds every year
%     Similar to Table 1 in Buehlmaier and Whited (2018) RFS, "Are financial constraints priced?"                     
%==========================================================================================================%

CHAR = xlsread('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\HedgeFundCharacters.xlsx','FundCharacteristics');
    %-- Creating dummies for Style
StyleD = dummyvar(CHAR(:,2));
StyleFF = StyleD(:,7);                                        % Fund of Fund

FundStyle = [CHAR(:,1) StyleFF];                              %[FundID Fund_of_Fund];
    %=== Hedge Fund Time Series Performance ===%
DATA_TS = xlsread('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\HedgeFundPerformance.xlsx');
    % DATA_TS = [FundID	Year Month Day Event Return	NAV	AUM]
    %-- Calculating Fund Flows --%
b = 1;
while b<=length(DATA_TS)
    n = length(find(DATA_TS(:,1)==DATA_TS(b,1)));
    DATA_TS(b,9) = -99;  %Indicator for missing value
    DATA_TS(b+1:b+n-1,9) = (DATA_TS(b+1:b+n-1,8) - (DATA_TS(b:b+n-2,8).*(1+DATA_TS(b+1:b+n-1,6))))./DATA_TS(b:b+n-2,8);
    b = b + n;
end
DATA_TS(DATA_TS(:,9)==-99,:) = []; % Delete missing values for Flow
        %DATA_TS = [FundID	Year Month Day Event Return	NAV	AUM Flow]
b=[]; n=[];
   %-- Cleaning data: Deleting funds less than 24 months --%
b=1; counter=0;
while b<=size(DATA_TS,1)
    counter = counter + 1;
    n = length(find(DATA_TS(:,1) == DATA_TS(b,1)));
    LifeSpan(b:b+n-1,1) = n;
    b = b + n;
end
DATA_TS(LifeSpan<=24,:) = [];                % Funds less than 24 months
DATA_TS(DATA_TS(:,9)>0.9,:)=[];              % Funds flow outliers >0.9
%-- Splitting funds into Non_FOF funds and FOF funds --%
for i = 1:size(DATA_TS,1)
    fof(i,1) = FundStyle(ismember(FundStyle(:,1),DATA_TS(i,1)),2);
end
DATA_IN = DATA_TS(fof==0,:);      % Non-FoF funds

%== Selecting both alive and liquidated hedge funds for every year starting from 1993 to 2021 ===%
T{1} = DATA_IN(DATA_IN(:,2)==1993,:);       % Those Non-FOF funds before 1993
Year = 1994:2021;
%W=vertcat([], T{1});
for t = 1:length(Year)
    T{t+1} = DATA_IN(DATA_IN(:,2)==Year(t),:);
end
%=== Compute Distance between every fund and the centriod for each cluster ===%
for k = 1:size(T,2)
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
    % X = [performance classfier];
    % %  GUI_YI X
    % X(:,1) =(X(:,1)-min(X(:,1)))/(max(X(:,1))-min(X(:,1)));
    % X(:,2) =(X(:,2)-min(X(:,2)))/(max(X(:,2))-min(X(:,2)));
    % X(:,3) =(X(:,3)-min(X(:,3)))/(max(X(:,3))-min(X(:,3)));
    % X(:,4) =(X(:,4)-min(X(:,4)))/(max(X(:,4))-min(X(:,4)));
    % if ~isnumeric(performance)
    %     error('Input data must be numeric.');
    % end
    % 
    % % Remove non-finite values
    % if any(~isfinite(performance(:)))
    %     error('Input data contains NaN or Inf values.');
    % end
    cleanedData = performance(all(isfinite(performance), 2), :);
    % Ensure data is real
    cleanedData = real(cleanedData);
    
    % Find rows with any NaN or Inf
    rowsWithNaNOrInf = any(~isfinite(performance), 2);
    
    % Display the indices of rows with NaN or Inf
    rowIndices = find(rowsWithNaNOrInf);
    
    classfier(rowIndices)=[];
    
    % Optional: Normalize data
    cleanperformance = normalize(cleanedData);
    
    X = [cleanperformance classfier];
    [~, ID_L, ID_M] = FcmGrnn_NonFoF_2SampleTtestA(X,ID,80);   
    % Measuring least Distressed nonFoFfunds performance
    RET_L(k,1) = mean(T{k}(ismember(T{k}(:,1),ID_L),6));
    STD_L(k,1) = std(T{k}(ismember(T{k}(:,1),ID_L),6));  
    AUM_L(k,1) = mean(T{k}(ismember(T{k}(:,1),ID_L),8));
    FLW_L(k,1) = mean(T{k}(ismember(T{k}(:,1),ID_L),9));
    % Measuring Most Distressed non-FOF funds performance
    RET_M(k,1) = mean(T{k}(ismember(T{k}(:,1),ID_M),6));
    STD_M(k,1) = std(T{k}(ismember(T{k}(:,1),ID_M),6));  
    AUM_M(k,1) = mean(T{k}(ismember(T{k}(:,1),ID_M),8));
    FLW_M(k,1) = mean(T{k}(ismember(T{k}(:,1),ID_M),9));
    % Two sample t-test
%     [H, Ret_P(k,1)] = ttest2(T{k}(ismember(T{k}(:,1),ID_L),6), T{k}(ismember(T{k}(:,1),ID_M),6));
%     [H, AuM_P(k,1)] = ttest2(T{k}(ismember(T{k}(:,1),ID_L),8), T{k}(ismember(T{k}(:,1),ID_M),8));
%     [H, Flw_P(k,1)] = ttest2(T{k}(ismember(T{k}(:,1),ID_L),9), T{k}(ismember(T{k}(:,1),ID_M),9)); 
    
    performance=[]; classfier=[]; ID=[]; ID_L=[]; ID_M=[]; X=[];
end
StatisticsCentriod = [(1993:2021)' RET_L RET_M  STD_L STD_M FLW_L FLW_M log(AUM_L) log(AUM_M)];

