function [T, G1, G2, Ig1, Ig2, Stat, Qn] = vwien400x(P, tmax, sw)
% ======================================================================
% Wiener process with a 3-stage gamma drift rate. (incomplete gamma)
% Uses Roger-style parameters: 
%    dX(t) = mu(t).dt + [s1(t) + s2].dB(t)   Ter variable.
% Parameters:
%    1.  mu: mean drift rate
%    2.  beta:  gamma dispersion
%    3.  m: number of gamma stages
%    4.  a: upper criterion
%    5.  z: starting point
%    6.  eta: drift standard deviation.
%    7.  sz:  starting point standard deviation.
%    8.  Ter: base time.
%    9.  st: Ter variability
%    10.  sigma1: variable component of variance.
%    11. sigma2: constant component of variance.
%    [T, G1, G2, Ig1, Ig2, Stat, Qn] = vwien400x(P, tmax, {switch})
%    Stat = [P1, P2, M1, M2, Skq1, Skq2]
%    sw: 0 or omitted, no Stat or Qn
%      : 1 Calculate Stat and Qn
%      : 2 Calculate Stat Qn and execution time.
% ====================================================================
NP = 11;
NS = 400; 
name = 'VWIEN400X: ';

if nargin < 3
   sw = 0;
   trace = 0;
end;
if sw == 2
   trace = 1;
   t0 = clock;
else
   trace = 0;
end;
if length(P) ~= NP
   disp([name, ' Incorrect length parameter vector, exiting...']);
   length(P)
   return;
end;

% Set me ! must agree with C-mex spacing.
% --------------------------------------      
h = tmax / NS;
epsilon = .01;
Stat = zeros(1,4);

% ------------------------------------------------------------------
% n-step probability mixture.
% Change the next two lines to change the quantisation or the range
% of the approximating truncated normal.
% ------------------------------------------------------------------
n_eta_step = 11;   % 19;  
n_sz_step =  11;   % 11;  
range = 4.5;       % 2.0;
U = ones(n_sz_step, 1);
    
delta = 2.0*range/n_eta_step;
Cbounds =[-10.0, -range + delta * [0:n_eta_step], 10];

% Generate normal mass between specified category bounds.
Pmass = diff(0.5 + 0.5 * erf(Cbounds/2^0.5));

% Rectangular mass
Rmass = U / n_sz_step ; 
 
% Pool mass outside plus/mivs range into end categories.
Pmass(2) = Pmass(2) + Pmass(1);
Pmass(n_eta_step+1) = Pmass(n_eta_step+1) + Pmass(n_eta_step+2);
Pmass = Pmass(2:n_eta_step+1);
Pstep = -range + delta * [0:n_eta_step-1]' + delta/2;
Rstep = [-(n_sz_step-1)/2:(n_sz_step-1)/2]' / (n_sz_step-1); 

% --------------------------------------------------------------------   
% Convert from s = .1 to sigma = 1.0 scaling.
% Parameters of wiensustn400x:
%      [v, beta, a, z, eta, sz, Ter, st, sigma1, sigma2]
% Scale (v, a, z, eta, sz, sigma1, sigma2) * 10 
% --------------------------------------------------------------------

Ix = [1,4,5,6,7,10,11];

P(Ix) = P(Ix) * 10.0; % Scale the sigma-dependent parameters.

P(Ix) = P(Ix) * 10.0; % Scale the sigma-dependent parameters.

mu = P(1);
beta = P(2);
m = P(3);
a = P(4);
z = P(5);
eta = P(6);
sz = P(7);
Ter = P(8);
st = P(9);
sigma1 = P(10);
sigma2 = P(11);

if (eta <= epsilon && sz <= epsilon)
    % ----------------
    % No variability.
    % ----------------
    a1 = a - z;    
    a2 = -z; 
elseif (eta > epsilon && sz <= epsilon)
    % -----------------------
    % Drift variability only.
    % -----------------------
    Mu = mu + Pstep * eta;
    a1 = a - z;   
    a2 = -z;
elseif (eta <= epsilon & sz > epsilon)    
    % ---------------------------
    % Criterion variability only.
    % ---------------------------
     Z = z + Rstep * sz;  
     A1 = (a * U - Z); 
     A2 = -Z;
     if (Z(1) < epsilon | Z(n_sz_step) > a - epsilon)
        disp([Z(1),Z(n_sz_step),a]);
        disp('Criterion range error');
        T=[]; G1=[]; G2=[]; Ig1=[]; Ig2=[]; Stat=[]; Qn=[];
        return;
     end;
else
    % ------------------------
    % Joint variability.
    % ------------------------
    Mu = (mu + Pstep * eta);  
    % Z =  z + Pstep * sz;<-- Normal
    Z = z + Rstep * sz;  % <- Rectangular
    A1 = (a * U - Z); 
    A2 = -Z;
    if (Z(1) < epsilon | Z(n_sz_step) > a - epsilon)
        disp([Z(1),Z(n_sz_step),P(2)]);
        disp('Criterion range error');
        T=[]; G1=[]; G2=[]; Ig1=[]; Ig2=[]; Stat=[]; Qn=[];
        return;
     end;   
end;    

if (eta <= epsilon && sz <= epsilon)
   % No variability to speak of.  
   [T, G1, G2, Ig1, Ig2, Stati] =  ...
        wiensustn400x([mu, beta, m, a1, a2, sigma1, sigma2], h);
   Pi = [mu, beta, a1, a2, sigma1, sigma2];
   Stat(1:2) = Stati(1:2); 
   Stat(3) = Stati(1) * Stati(3);
   Stat(4) = Stati(2) * Stati(4);
   %plot(T, G1, T, G2)
   %Stati
   %pause
else
   % ------------------------------------------------------------------
   % Some variability, do mixtures.
   % ------------------------------------------------------------------
    G1 = zeros(1,NS);
    G2 = zeros(1,NS);    
    Ig1 = zeros(1,NS);    
    Ig2 = zeros(1,NS);
    T = zeros(1, NS);
    Stat1 = 0;
    Stat2 = 0;
    Stat3 = 0;
    Stat4 = 0;
    Stat = zeros(1, 4);
    if (eta > epsilon && sz <= epsilon)
         % ----------------------------------------
         % Variability in drift rates.
         % ----------------------------------------
        %parfor i = 1:n_eta_step
        for i = 1:n_eta_step
   	     [Ti, Gi1, Gi2, Igi1, Igi2, Stati] =  ...
        	    wiensustn400x([Mu(i), beta, m, a1, a2, sigma1, sigma2], h);
             T = T + Ti; 
             G1 = G1 + Pmass(i) * Gi1;
             G2 = G2 + Pmass(i) * Gi2;
             Ig1 = Ig1 + Pmass(i) * Igi1;
             Ig2 = Ig2 + Pmass(i) * Igi2;
             Stat1 = Stat1 + Pmass(i) * Stati(1);
             Stat2 = Stat2 + Pmass(i) * Stati(2);            
             Stat3 = Stat3 + Pmass(i) * Stati(1) * Stati(3);
             Stat4 = Stat4 + Pmass(i) * Stati(2) * Stati(4);
             %Stati
             %plot(T, Gi1, T, Gi2)
             %pause
         end;
         Stat = [Stat1, Stat2, Stat3, Stat4];
         T = T / n_eta_step; % Because of parfor
         %disp('Calling 2')
    elseif (eta <= epsilon & sz > epsilon)
         for i = 1:n_sz_step
             % [A1(i), A2(i)]
	     [T, Gi1, Gi2, Igi1, Igi2, Stati] =  ...
        	 wiensustn400x([mu, beta, m, A1(i), A2(i), sigma1, sigma2], h);
             G1 = G1 + Rmass(i) * Gi1;
             G2 = G2 + Rmass(i) * Gi2;
             Ig1 = Ig1 + Rmass(i) * Igi1;
             Ig2 = Ig2 + Rmass(i) * Igi2;
             Stat(1) = Stat(1) + Rmass(i) * Stati(1);
             Stat(2) = Stat(2) + Rmass(i) * Stati(2);
             Stat(3) = Stat(3) + Rmass(i) * Stati(1) * Stati(3);
             Stat(4) = Stat(4) + Rmass(i) * Stati(2) * Stati(4);
         end;
         %disp('Calling 3')
    else
         % -------------------------------------------------------
         % Joint variability in both drift and starting position.
         % -------------------------------------------------------

        for i = 1:n_eta_step
             for j = 1:n_sz_step
		[T, Gij1, Gij2, Igij1, Igij2, Statij] =  ...
                     wiensustn400x([Mu(i), beta, m, A1(i), A2(i), sigma1, sigma2], h);
                 G1 = G1 + Pmass(i) * Rmass(j) * Gij1;
                 G2 = G2 + Pmass(i) * Rmass(j) * Gij2;
                 Ig1 = Ig1 + Pmass(i) * Rmass(j) * Igij1;
                 Ig2 = Ig2 + Pmass(i) * Rmass(j) * Igij2;
                 Stat(1) = Stat(1) + Pmass(i) * Rmass(j) * Statij(1);
                 Stat(2) = Stat(2) + Pmass(i) * Rmass(j) * Statij(2);
                 Stat(3) = Stat(3) + Pmass(i) * Rmass(j) * Statij(1) * Statij(3);
                 Stat(4) = Stat(4) + Pmass(i) * Rmass(j) * Statij(2) * Statij(4);
             end;
         end;
         %disp('Calling 4');
    end;
end;
Stat(3) = Stat(3) / max(Stat(1), epsilon) + Ter;    % E[T1]
Stat(4) = Stat(4) / max(Stat(2), epsilon) + Ter;    % E[T2]
T = T + Ter;

% --------------------------------------------------------
% Convolve with Ter.
% --------------------------------------------------------
if st > 2 * h
   m = round(st/h);
   n = length(T);
   fe = ones(1, m) / m;
   Ig1 = conv(Ig1, fe); 
   Ig2 = conv(Ig2, fe);
   Ig1 = Ig1(1:n);
   Ig2 = Ig2(1:n);
end;

% --------- mods --------------
m1 = max(Ig1);
m2 = max(Ig2);
Ig1 = Ig1 / (m1 + m2);
Ig2 = Ig2 / (m1 + m2);
% -----------------------------
% Calculate quantiles.
if sw
   Qf5 = [.1, .3, .5, .7, .9]; % Summary quantiles (can be changed).

   % Normalise distributions, select subrange spanning quantiles.
   F1 = Ig1 / max(Ig1);
   F2 = Ig2 / max(Ig2);
   Dummy5 = [.001,.002,.003,.004,.005];

   I1 = (F1 >= .04 & F1 <= .96); 
   if min(diff(F1(I1))) <= 0 
        Q1 = Dummy5;
        disp('Cannot compute F1 quantiles.');    
   else
        Q1=interp1(F1(I1)', T(I1)', Qf5);
   end;

   I2 = (F2 >= .04 & F2 <= .96);  
   if min(diff(F2(I2))) <= 0 
        Q2 = Dummy5;
        disp('Cannot compute F2 quantiles.');    
   else
        Q2=interp1(F2(I2)', T(I2)', Qf5);
   end;
   Qn = [Q1',Q2'];
else
   Qn = [];
end;

if sw == 2
     fprintf('ET = %10.3f \n',etime(clock, t0));
end;

