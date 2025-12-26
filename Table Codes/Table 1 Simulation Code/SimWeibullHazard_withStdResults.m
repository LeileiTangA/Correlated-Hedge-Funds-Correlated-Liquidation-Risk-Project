tic
rng(12); 
M = 1000;           % # of simulations
matrixSize = 800;  % # of hedge funds
Beta = [1 -1]';
density = 0.1:0.2:0.9;      % Sparse density
for i = 1:M
    [Coef_Corr(i,:), ~, Coef_Inde(i,:),Coef_CorrZ(i,:), Coef_CorrZPD(i,:) ] = SimCorrWeibullCox_WithStd(matrixSize, Beta, density(2));
end
Result = [mean(Coef_Corr)  mean(Coef_Inde)            mean(Coef_CorrZ)                 mean(Coef_CorrZPD)]
%         AllCorrleated    AllIndependent    AllCorrelate2distressed_fund    AllIndependt_but_correlate2distressed_fund
toc
Result_Std = [std(Coef_Corr)  std(Coef_Inde)            std(Coef_CorrZ)                 std(Coef_CorrZPD)]

% Result =
% 
%   Columns 1 through 9
% 
%     2.0924   -2.0873    1.0082   -1.0093    0.8941   -0.8897    1.3351    1.0067   -1.0068
% 
%   Column 10
% 
%     1.5118
% 
% Elapsed time is 2269.499437 seconds.
% 
% Result_Std =
% 
%   Columns 1 through 9
% 
%     0.6710    0.6706    0.0647    0.0860    0.0899    0.0999    0.2063    0.0588    0.0713
% 
%   Column 10
% 
%     0.1349
function [Coef_Corr, Coef_Spar,Coef_Inde,Coef_CorrZ, Coef_CorrZPD]= SimCorrWeibullCox_WithStd(matrixSize, Beta, density)
% This function is to simulated correlated survival time for Weibull-Cox hazard regression
% Input: matrixSize =  #_of_hedgefunds; Beta = True Coefficients like Beta = [1 1 -1]';
% Output: Coef = Estimated Coeffcients

%-- Generate corrleated survival time from Weibull distribution--%
%-- Generate a random symmetrical matrix with a diagonal line of 1s --%
% For any square matrix A, A' * A is positive semi-definite, and rank(A' * A) is equal to rank(A) . 
% Matrices are invertible if they have full rank. So all we have to do is generate an initial random matrix with full rank 
% and we can then easily find a positive semi-definite matrix derived from it. Nearly all random matrices are full rank, 
% so the loop I show will almost always only iterate once and is very very unlikely to need more than a very small number of iterations.

%matrixSize = 1500;
while true
  A = rand(matrixSize, matrixSize);
  if rank(A) == matrixSize; break; end    %will be true nearly all the time
end

A = corrcov(A' * A);
%== Simulate correlated Survival Time ==%
[V,D] = eig(A);
V12 = V*sqrt(D)*V';
Yrnd = normrnd(0,1,[matrixSize,1]);
Y = V12*Yrnd;
U = normcdf(Y);
%-- Define Explanatory variables and Coefficients --%
%X = [ones(matrixSize,1) normrnd(0,1,[matrixSize 1]) chi2rnd(1,matrixSize,1)];
X = [normrnd(0,1,[matrixSize 1]) chi2rnd(1,matrixSize,1)];
    % X = [normalDist Chi2Dist(with 1 degree of freedom)]
%Beta = [1 1 -1]';  
Lamda = exp(X*Beta);
shape = 2; scale = 0.1;% Relevant parameters for Weibull Distribution
% Generate latent randomly correlated survival time
T_Surv = ((- log(U))./(scale * Lamda)).^(1 / shape);
% Gnerate Independent corrleated survival time
T_Inde = ((- log(normcdf(Yrnd)))./(scale * Lamda)).^(1 / shape);
% Generate latent Sparsely randomly corrleated survial time
B = generatesparseSPDmatrix(matrixSize,density);
[VB,DB] = eigs(B);
V12B = VB*sqrt(DB)*VB';
YB = V12B*Yrnd;
UB = normcdf(YB);
T_Spar = ((- log(UB))./(scale * Lamda)).^(1 / shape);

%-- Define censoring time and event indicator --%
rateC = 4;  % A fund is expected to last an average of Four years
T_Cens = exprnd(rateC, matrixSize,1); % Censoring times
% Correlated survival time
Time = min(T_Surv,T_Cens);
Event = (T_Surv <= T_Cens);
% Independent survival time
Time_Inde = min(T_Inde,T_Cens);
Event_Inde = (T_Inde <= T_Cens);
% Sparse corrleated survival time
Time_Spar = min(T_Spar,T_Cens);
Event_Spar = (T_Spar <= T_Cens);

%== Fit for Weibull Hazard regression ==%
warning('off')
if rank(X'*X) ~= size(X,2); return; end
model = wc_train(X, Time, Event);
Coef_Corr = model.beta';

modelInde = wc_train(X, Time_Inde, Event_Inde);
Coef_Inde = modelInde.beta';

modelSpar = wc_train(X, Time_Spar, Event_Spar);
Coef_Spar = modelSpar.beta';

%=== The following is to conduct the same above. However, the simulated corrleation matrix is conditional on another fund, saying 'Fund A'. 
% See if conditional correlation matrix (Partical corrleation matrix) can be improve the survival estimation 
Z = rand(matrixSize+1, matrixSize+1); % This '+1' is 'Fund A'
% Calcuating partial corrleation matrix conditional on the last one
PZ = partialcorr(Z(1:end-1,1:end-1),Z(1:end-1,end));
RHO = Z(1:end-1,end);   % This is correlations between each fund with the last one. It will be used as an explanatory variable in the hazard model
PZA = FullRankCorr(PZ); % PZA is positive & definite
%== Simulate correlated Survival Time ==%
[VZ,DZ] = eig(PZA);
VZ12 = VZ*sqrt(DZ)*VZ';
YZ = real(VZ12)*Yrnd;
UZ = normcdf(YZ);
%-- Define Explanatory variables and Coefficients --%
XZ = [X RHO];
    % X = [normalDist Chi2Dist(with 1 degree of freedom) RHO]
LamdaZ = exp(XZ*[Beta; 1.5]);
% Generate latent randomly correlated survival time
T_SurvZ = ((- log(UZ))./(scale * LamdaZ)).^(1 / shape);
TimeZ = min(T_SurvZ,T_Cens);
EventZ = (T_SurvZ <= T_Cens);
if rank(XZ'*XZ) ~= size(XZ,2); return; end
model = wc_train(XZ, TimeZ, EventZ);
Coef_CorrZ = model.beta';

%=== The following is to conduct the same above. However, all funds are indepdendent excepf for fund A. 
% See if conditional correlation matrix (Partical corrleation matrix) can be improve the survival estimation 
Zrho = rand(matrixSize,1); % This is correlation for all funds with 'Fund A'
ZA = eye(matrixSize+1);    % Correlation matrix for all funds + 'Fund A'
ZA(matrixSize+1,1:matrixSize)=Zrho;
ZA(1:matrixSize,matrixSize+1)=Zrho; % Correlation matrix for all funds + 'Fund A'
ZAPD = FullRankCorr(ZA);  % PZAPD is positive & definite
% Calcuating partial corrleation matrix conditional on the last one
PZAPD = partialcorr(ZAPD(1:end-1,1:end-1),ZAPD(1:end-1,end));
%== Simulate correlated Survival Time ==%
[VZPD,DZPD] = eig(PZAPD);
VZ12PD = VZPD*sqrt(DZPD)*VZPD';
YZPD = real(VZ12PD)*Yrnd;
UZPD = normcdf(YZPD);
%-- Define Explanatory variables and Coefficients --%
XZPD = [X Zrho];
    % X = [normalDist Chi2Dist(with 1 degree of freedom) RHO] 
LamdaZPD = exp(XZPD*[Beta; 1.5]);
% Generate latent randomly correlated survival time
T_SurvZPD = ((- log(UZPD))./(scale * LamdaZPD)).^(1 / shape);
TimeZPD = min(T_SurvZPD,T_Cens);
EventZPD = (T_SurvZPD <= T_Cens);
if rank(XZPD'*XZPD) ~= size(XZPD,2); return; end
modelPD = wc_train(XZPD, TimeZPD, EventZPD);
Coef_CorrZPD = modelPD.beta';
end

function A_PD = FullRankCorr(A)
% https://uk.mathworks.com/matlabcentral/answers/320134-make-sample-covariance-correlation-matrix-positive-definite
% Sample covariance and correlation matrices are by definition positive semi-definite (PSD), not PD. 
% Semi-positive definiteness occurs because you have some eigenvalues of your matrix being zero (positive definiteness guarantees all your eigenvalues are positive).
% If you correlation matrix is not PD ("p" does not equal to zero) means that most probably have collinearities between the columns of your correlation matrix, 
% those collinearities materializing in zero eigenvalues and causing issues with any functions that expect a PD matrix.
 
% To fix this the easiest way will be to do calculate the eigen-decomposition of your matrix and set the "problematic/close to zero" eigenvalues 
% to a fixed non-zero "small" value. That can be easily achieved by the following code, given your initial correlation matrix "A":
[V,D] = eig(A);       % Calculate the eigendecomposition of your matrix (A = V*D*V') 
                        % where "D" is a diagonal matrix holding the eigenvalues of your matrix "A"
d= diag(D);           % Get the eigenvalues in a vector "d" 
d(d <= 0.1) = 0.1;  % Set any eigenvalues that are lower than threshold "TH" ("TH" here being 
                        % equal to 1e-7) to a fixed non-zero "small" value (here assumed equal to 1e-7)
D_c = diag(d);        % Built the "corrected" diagonal matrix "D_c"
A_PD = V*D_c*V';      % Recalculate your matrix "A" in its PD variant "A_PD"
end

function A = generatesparseSPDmatrix(n,density)
% Generate a sparse n x n symmetric, positive definite matrix with approximately density*n*n non zeros
A = sprandsym(n,density); % generate a random n x n matrix
% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix is symmetric positive definite, which can be ensured by adding I
A = A + speye(n);
A = corrcov(A'*A);
end


