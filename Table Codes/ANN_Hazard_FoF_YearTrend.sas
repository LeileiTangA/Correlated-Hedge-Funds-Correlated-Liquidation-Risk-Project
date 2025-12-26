PROC IMPORT OUT=ANNHazard 
            DATAFILE= "C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\YearTrend\SAS4AML_YearTrend.xlsx" 
            DBMS=xlsx REPLACE; GETNAMES=YES;
RUN;
/* Survival analysis for FCM-ANN hedge fund liquidation: June 11 2020 */
data SurviveAnalysis; 
	set ANNHazard;		
	weight=1;
	/* == Save data for aML== */
	file 'C:\Users\ias05106\OneDrive - University of Strathclyde\Desktop\Updated Hedge Fund Project\YearTrend\ANNHazard_YearTrend.raw'; 	/* Liquidation hazard */
	put FundID  Numint					        /*Control Variables*/
		weight						            /*Level 1 Variable*/
		Event DropID Year1992 Lower Upper; 				/*Level 2 Variable*/
								
		array t(*) t1-t30;
		array DistL(*) DistL1-DistL30;
		array DistH(*) DistH1-DistH30;
		array Ret(*)   Ret1-Ret30;
		array Volt(*)  Volt1-Volt30;
		array Aum(*)   Aum1-Aum30;
		array Flow(*)  Flow1-Flow30;
		array VIX(*)   VIX1-VIX30;
		array CSprd(*) CSprd1-CSprd30;
		array AggLq(*) AggLq1-AggLq30;
		array InnLq(*) InnLq1-InnLq30;
		array TrdLq(*) TrdLq1-TrdLq30;
		array CRSP(*)  CRSP1-CRSP30;
		DO i=1 to Numint;
		put t(i) DistL(i) DistH(i) Ret(i) Volt(i) Aum(i) Flow(i)    /*Level 3 Variable*/
		         VIX(i) CSprd(i) AggLq(i) InnLq(i) TrdLq(i) CRSP(i);
	    END;
run;
