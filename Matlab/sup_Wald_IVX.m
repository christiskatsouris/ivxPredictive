%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R Script Details:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script name: supWaldIVX_single_predictor.R

% Program aim: This R program implements the sup Wald IVX test for structural break detection 
% on the slopes of a linear predictive regression model with two predictors.  

% written by: 

% Christis G. Katsouris (August 2021)
% Department of Economics
% University of Southampton
% Southampton, United Kingdom


function [sup_Wald_IVX] = estimate_sup_Wald_IVX_quantile( ysim, xsim );

outputVector = [];

% Inputs:
% Simulated DGP under the null hypothesis 
% Model: no intercept and 2 predictors

% BEGIN OF FUNCTION 
% trimming coefficient
pi0 = 0.15;
p = 2;
n = length(ysim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the simulated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlag = xsim(1:n-1,:);
xt   = xsim(2:n,:);
y    = ysim(2:n,:);
nn   = length(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit the LUR specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rn = zeros(p,p);
for i = 1:p
   rn(i,i) = regress( xt(:,i),xlag(:,i) );
end

% autoregressive residual estimation
u = xt - xlag*rn;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the IVX instrument  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
K = 1;
n  = nn-K+1;
Rz = ( 1 - 1/(nn^0.95) )*eye(p); 
diffx = xt-xlag; 
z = zeros(nn,p);
z(1,:) = diffx(1,:);
    
for i=2:nn
  z(i,:) = z(i-1,:)*Rz + diffx(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop for the supremum functional  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = length( y ); 
lower = t*pi0;
upper = t*(1-pi0);

% qs    = sort( v2(1:t-1, 1) ,1 );
% qss   = qs( (lower+1):upper, : )
% dim   = length( qss );

lower_bound = round(lower);
upper_bound = round(upper);
sequence    = ( lower_bound : upper_bound );
dim_seq     = ( upper_bound - lower_bound );

% Wald IVX vector that stores the sequence of statistics
% for the nonstationary quantile predictive regression model

Wald_IVX_vector = zeros(dim_seq ,1);
s = 1;

while ( s <= dim_seq  );
   
    k = sequence( s )
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: Obtain the estimate of the variance of OLS regression 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
    xlag1 = xlag;
    xlag1( (k+1):nn, : ) = 0;
    
    xlag2 = xlag;
    xlag2(  1:k, : ) = 0;
    
    X_matrix = [ xlag1 xlag2 ];
        
    [Aols,bhat,epshat] = regress( y, X_matrix ); 
    covepshat = ( epshat'*epshat ) / nn;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2: Obtain the estimates of the IVX estimators
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % We define with Z the IVX instrument    
    Z = [ zeros(1,p) ; z(1:n-1,:) ];
    
    % Estimate the Z1 (before break) and Z2 (after break) matrices     
    Z1 = Z;
    Z1( (k+1):nn, : ) = 0;
    
    Z2 = Z;
    Z2( 1:k, : ) = 0;
    
    beta1_ivx = y'*Z1*pinv(xlag1'*Z1);
    beta2_ivx = y'*Z2*pinv(xlag2'*Z2);
    beta_ivx_distance = ( beta1_ivx - beta2_ivx );
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimation of Q1 and Q2 matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Z1_mean   = mean( Z1( 1:k, : ) ); 
    M1_matrix = Z1'*Z1*covepshat - k*Z1_mean'*Z1_mean*FM;
    Q1_matrix = pinv(Z1'*xlag1)*( M1_matrix )*pinv(xlag1'*Z1);
   
    Z2_mean   = mean( Z2( (k+1):nn, : ) ); 
    M2_matrix = Z2'*Z2*covepshat - ( n - k )*Z2_mean'*Z2_mean*FM;
    Q2_matrix = pinv(Z2'*xlag2)*( M2_matrix )*pinv(xlag2'*Z2);
    
    Q_matrix       = ( Q1_matrix + Q2_matrix );
    Wald_IVX_break = beta_ivx_distance*pinv( Q_matrix )*(beta_ivx_distance'); 
    
    Wald_IVX_vector( s, 1 ) = Wald_IVX_break;
    
% end-of estimation of sequence of Wald-IVX statistics     
s = s + 1;
end

% Obtain the sup-Wald statistic
sup_Wald_IVX = max(Wald_IVX_vector);

end
%END-OF-FUNCTION
