%==== Liquidated Funds  ====%
Performance_Liq = readtable('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\Lipper TASS 2021\tassfarc3\ProductPerformance.txt','headerlines',0,'delimiter',',','ReadVariableNames', true);
          % Sorting date in order to get Year Month infomation
t = string(Performance_Liq.Date);
d = datetime(t) ;
d.Format = 'yyyy-MM-dd';
column1 = cellstr(d);         % Cell array of strings.
% d.Format = 'hh:mm:ss';
% column2 = cellstr(d);
[Year,Month,Day]=ymd(datetime(column1));
          % Get data ready
Event_Liq = zeros(size(Performance_Liq,1),1);   % Fund Live=Liquidation indicator 0 = live, 1 = Liquidation
          % Assigning the last month for each fund event indicator to be 1,i.e., fund had been liquidated
b=1;
counter = 0;
while b <= size(Performance_Liq,1)
    counter = counter + 1;
    n = length(find(Performance_Liq.ProductReference == Performance_Liq.ProductReference(b,1)));
    Event_Liq(b+n-1,1) = 1;
    b = b + n;
end
% PerformanceMeasure_Liq = [Performance_Liq.ProductReference Year Month Day Event_Liq...
%                            Performance_Liq.RateOfReturn/100 Performance_Liq.NAV Performance_Liq.EstimatedAssets];
% Header={'FundID','Year','Month','Day','Event', 'Return' , 'NAV', 'AUM'};
% PM = [PerformanceMeasure_Live; PerformanceMeasure_Liq];

PM = [Performance_Liq.ProductReference Year Month Day Event_Liq];                   
CHAR = xlsread('C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\HedgeFundCharacters.xlsx','FundCharacteristics');
     %-- Creating dummies for Style
StyleD = dummyvar(CHAR(:,2));
StyleFF = StyleD(:,7);     % Fund-of-Fund Dummy
FoF = [CHAR(:,1) StyleFF]; %[FundID FoF_dummy]
PM(:,6)=zeros(size(PM,1),1);
for i = 1:size(PM,1)
    f=find(FoF(:,1)==PM(i,1));
    if length(f)>0
        PM(i,6)=FoF(f,2);
    else
        PM(i,6)=0;
    end
    f=[];
end
% Counting the number of liquidated Fund-of-Funds and other style hedge funds each year from 1994 to 2018.
Year = 1993:2021;
for i = 1:length(Year)
    Liquid_FoF(i,1) = sum(PM(PM(:,2)==Year(i) & PM(:,6)==1 & PM(:,5)==1,5));% Fund of Funds
    Liquid_Oth(i,1) = sum(PM(PM(:,2)==Year(i) & PM(:,6)==0 & PM(:,5)==1,5));% Other types of funds
end
Q = [Year' Liquid_FoF Liquid_Oth];
figure
hold on
set(gca,'YAxisLocation','origin')
for i=1:29
    plot([0,log(Q(i,2))],[Q(i,1)-0.5,Q(i,1)-0.5],'g-','linewidth',2);
    text(log(Q(i,2))+0.5,1990.5+i,num2str(Q(i,2)),'vert','bottom','horiz','center'); 
    plot([-log(Q(i,3)),0],[Q(i,1)-0.5,Q(i,1)-0.5],'r-','linewidth',2);
    text(-log(Q(i,3))-0.5,1990.5+i,num2str(Q(i,3)),'vert','bottom','horiz','center'); 
end
%set(gca,'XtickLabel',{'(1200)','(800)','(400)','(100)','0','50','100' '200','400'})
set(gca,'XtickLabel',{[]})
set(gca,'Ytick',[1993:2021])
yticklabels({'93','94','95','96','97','98','99','00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18' ,'19','20','21'})
set(gca,'TickLength',[0;0])
legend('FoFs','Hedge Funds')
title('Figure 1: Patterns of Liquidated FoFs(right) vs. Hedge Funds(left)')
xlabel('Annual number of liquidated funds')
%ylabel('Year')

    












