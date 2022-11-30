%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date of modification: 15 November 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = 100;
M = 100;
n = 200;
beta1 = 0.25;
beta2 = 0.75;

pi0 = 0.15;
p = 2;

gamma = 1;
rho   = 0.3;

%[ ysim , xsim ]  = simulate_null( n, beta1, beta2 );

% Given values of beta1 and beta2 (no model intercept)
[ysim, xsim] = simulate_null1( n , beta1, beta2, c1, c2, gamma, rho );

ysim 
xsim

xlag = xsim(1:n-1,:);
xt   = xsim(2:n,:);
y    = ysim(2:n,:);
nn   = length(y);

Wald_IVX          = estimate_sup_Wald_IVX( ysim, xsim );
Wald_IVX_quantile = estimate_sup_Wald_IVX_quantile( ysim, xsim );

% Fixed regressor bootstrap (need to modify this step)
yt_star = zeros(n-1,B);

sup_Wald_IVX_matrix           = zeros(M,1);
sup_Wald_IVX_matrix_bootstrap = zeros(M,1);

empirical_size = 0;
p_value_boot   = 0;
count = 0;
ysim  = 0;
xsim  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MC simulation step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;

for i = 1:M
%begin-of-MC-loop
   
   [ ysim , xsim ]  = simulate_null1( n , beta1, beta2, c1, c2, gamma, rho );
   xlag = xsim(1:n-1,:);
   xt   = xsim(2:n,:);
   y    = ysim(2:n,:);
   nn   = length(y);

   Wald_IVX = estimate_sup_Wald_IVX( ysim, xsim );
   sup_Wald_IVX_matrix( i, 1 ) = Wald_IVX;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Bootstrap step (Fixed regressor bootstrap)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   yt_star = zeros( nn,B );

   for b = 1:B
      [Aols,bhat,epshat] = regress( y, xlag ); 
      kappa      = normrnd(0,1,[nn,1]);
      u_hat_star = ( epshat.*kappa); 
      Aols(1,1)*xlag( :,1) + Aols(2,1)*xlag( :,2) + u_hat_star;  
   end   

   bootstrap_Wald_IVX_matrix = zeros( B , 1 );

   R = 8;
   parfor (j = 1:B, R)
       ysim = yt_star( : , j ) ;
       estimation_sup_Wald_IVX_step     = estimate_sup_Wald_IVX( ysim, xsim );
       bootstrap_Wald_IVX_matrix( j, 1) = estimation_sup_Wald_IVX_step;
   end

   count_boot     = 0;
   for s = 1:B
       if  ( bootstrap_Wald_IVX_matrix( s, 1) > Wald_IVX );
           count_boot = count_boot + 1;
       end
   end    

   p_value_boot =  count_boot / B;

   %the_boot_statistic = sort( bootstrap_Wald_IVX_matrix );
   %quantile           = round(0.95*B);
   %boot_IVX_statistic = the_boot_statistic( quantile,  1 );

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   sup_Wald_IVX_matrix_bootstrap( i, 1 ) = p_value_boot;
   ysim = 0;
   xsim = 0;

%end-of-MC-loop
