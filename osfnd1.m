function [G2, AIC, BIC, Pred] = osfnd1(Pvar, Pfix, Sel, Dmatrix, trace)
% ========================================================================
% Fit standard diffusion model to the OSF19 data using G2.
% Drift bias version with drift criterion c
% Data matrix "Dmatrix" must be 12 x 14 matrix. Uses median error RT
% if quantiles can't be computed.
% Usage:
%    [G2, AIC, BIC, Pred] = osfnd1(Pvar,Pfix, Sel, Dmatrix, {trace});
%    
%    P = [as, aa,  vh, ve,  c,  pz, eta, szs, sza, Ter  st]
%          1   2   3    4   5   6   7   8     9    10   11
%  "trace" is optional parameter; if 1, gives information about 
%  parameters that are out-of-bounds
%  Starting point parameterized as a proportion -1:1 of distance from 
% midpoint to boundary
% =========================================================================

Q1 = Dmatrix(:, [3,10]); % Pull out .1 quantile
ValidQ1 = Q1(Dmatrix(:, [3,10]) > 0); % Valid Q1 for min criterion.
AllValid = [all(Dmatrix(:,3:7) > 0, 2), all(Dmatrix(:,10:14) > 0, 2)];
AnyValid = [any(Dmatrix(:,3:7) > 0, 2), any(Dmatrix(:,10:14) > 0, 2)];
t0 = clock;
name = 'osfnd1: ';
errmg1 = 'Incorrect number of parameters for model, exiting...';
errmg2 = 'Incorrect length selector vector, exiting...';
errmg3 = 'Requires (12,14) size dat matrix, exiting...';
if nargin < 5
    trace = 0;
end;

% -----------------------------------
% Manual setting of tmax needed here.
% -----------------------------------
tmax = 3.0;
np = 11; 
nsteps = 300; % Number of time steps.
ncond = 12;

lp = length(Pvar) + length(Pfix);
if lp ~= np
     [name, errmg1], length(Pvar), length(Pfix), return;
end;
if length(Sel) ~= np
     [name, errmg2], length(Sel), return;
end;    
   
% Assemble parameter vector.
P = zeros(1,np);
P(Sel==1) = Pvar;
P(Sel==0) = Pfix;

Psave = P;
save Psave Psave
if trace
   P
end;
if any(size(Dmatrix) ~= [ncond,14])
     [name, errmg3], size(Dmatrix), return;
else
   Pf1 = [0,.1,.3,.5,.7,.9,1.0];
   Pred = zeros(ncond, 14);  % Predicted values.
end;

%N = 2800 / ncond; % Approximately
Ni = sum(Dmatrix(:, [2,9]), 2);
N = sum(Ni);

% Parameter bounds (absolute)
% Constrain Ter to be smaller than the smallest first data quantile.
epsn = .001;
as = P(1);
aa = P(2);
vh = P(3);
ve = P(4);
c = P(5);
pz = P(6);
eta = P(7);
szs = P(8);
sza = P(9);
Ter = P(10);
st = P(11);

% Hi F Hard  (Speed)
% Hi F Easy
% Eq Hard
% Eq Easy
% Lo F Hard
% Lo F Easy
% Repeat (Accuracy)

% Parameters for 12 conditions
Pi = [as, vh+c, (1+pz)*as/2, eta, szs, Ter, st;
      as, ve+c, (1+pz)*as/2, eta, szs, Ter, st;
      as, vh ,  as/2, eta, szs, Ter, st;
      as, ve,   as/2, eta, szs, Ter, st;
      as, vh-c, (1-pz)*as/2, eta, szs, Ter, st;
      as, ve-c, (1-pz)*as/2, eta, szs, Ter, st;
      aa, vh+c,  (1+pz)*aa/2, eta, sza, Ter, st;
      aa, ve+c, (1+pz)*aa/2, eta, sza, Ter, st;
      aa, vh,   aa/2, eta, sza, Ter, st;
      aa, ve,   aa/2, eta, sza, Ter, st;
      aa, vh-c, (1-pz)*aa/2, eta, sza, Ter, st;
      aa, ve-c, (1-pz)*aa/2, eta, sza, Ter, st];

q1min = min(ValidQ1);

%Zs = as - [(1-pzh)*as/2, (1-pze)*as/2, (1-pzl)*as/2];
% Minimum distance of any starting value to either boundary.
ZBounds = [(1-pz)*as/2, (1+pz)*as/2];
ZBounda = [(1-pz)*aa/2, (1+pz)*aa/2];

% Minimum distance of any starting value to either boundary.
szmaxs = min(ZBounds) - 0.001;
szmaxa = min(ZBounda) - 0.001;
mina = min(as, aa);

% -------------------------------------------------------------------------------------
%      as   aa    vh      ve    c      pzh     eta      szs   sza     Te      st ---------------------------------------------------------------------------------------
Ub =  [.28, .28,  0.50,  0.50,  0.5,  1.0,    0.5,  szmaxs,     szmaxa,   q1min,  0.5];
Lb =  [.05, .05, -0.50  -0.50, -0.5, -1.0,      0,      0,        0           0.1, 0 ];
Pub = [.25, .25,  0.45,  0.45,  0.45,  0.95,  0.45, szmaxs-epsn,szmaxs-epsn, q1min-0.05, .45];
Plb = [.06, .06, -0.45, -0.45, -0.25, -0.95     0,      0,        0          0.15, 0];

if trace
    P, Ub, Lb
    max(P - Ub, 0)
    max(Lb - P, 0)
end    

if any(P - Ub > 0) || any(Lb - P > 0)
   G2 = 1e5 + ...
         1e3 * (sum(max(P - Pub, 0).^2) + sum(max(Plb - P).^2));
   AIC = 1e8;
   BIC = 1e8;     
else
   % ------------------------
   % Set me!
   % ------------------------
   h = tmax/nsteps;
                                                     % Adjust for step size.  (Set h here!)
   Qni = [3:7,10:14];                                % Index the quantiles in data matrix.
   Sni = [1,2,8,9];                                  % Index of summary stats in data matrix.

   G2 = 0;
   ObsBuff = zeros(5,12);
   PredBuff = zeros(5,12);
   for i = 1:ncond
        %i
        %Pi(i,:)
        [t,g1,g2,ig1,ig2,StPr,QnPr] = etawiener(Pi(i,:), h, tmax, 1);
        %StPr
        %QnPr
        %plot(t, g1, t, g2)
        %pause
        Pred(i,Qni) = [QnPr(:,1)',QnPr(:,2)'];  % transpose quantiles to Pred matrix.
        Pred(i,Sni) = StPr([1,3,2,4]);
        % Predicted distribution functions linearly interpolated at the quantiles.
        % Trap quantile errors caused by Ter overrun. 
        Qnrlo = floor((Dmatrix(i,Qni) - Ter)/ h) + 1;
        Qnrlo = max(Qnrlo,1);
        Tlo = [t(Qnrlo(1:5)),t(Qnrlo(6:10))];
        Tfrac = (Dmatrix(i,Qni) - Tlo) ./ h;
        Qnrhi = ceil((Dmatrix(i,Qni) - Ter) / h) + 1;

        % Predicted values weighted implicitly by probability.
        Iglo = [ig1(Qnrlo(1:5)),    ig2(Qnrlo(6:10))];
        Ighi = [ig1(Qnrlo(1:5)+1),  ig2(Qnrlo(6:10)+1)];
        Ig = [Iglo(1:5) + Tfrac(1:5) .* (Ighi(1:5) - Iglo(1:5)), ...
              Iglo(6:10) + Tfrac(6:10) .* (Ighi(6:10) - Iglo(6:10))];    % (10 x 1)
        % Objective function.
        % Modified version to make a proper chi-square.
        %i
        %Dmatrix(i, [1,8])
        %StPr 
        PredProb1 = diff([0, Ig(1:5), max(ig1(length(t)), StPr(1))]);
        ObsProb1 = diff(Dmatrix(i,1) * Pf1);
        PredProb2 = diff([0, Ig(6:10), max(ig2(length(t)), StPr(2))]);
        ObsProb2 =  diff(Dmatrix(i,8) * Pf1);
        if AllValid(i,2)   % No missing quantiles
            PredProb = [PredProb1, PredProb2];
            ObsProb = [ObsProb1, ObsProb2];
        elseif AnyValid(i,2) % Use median
            PredProb2Pool = sum(PredProb2(4:6));
            ObsProb2Pool = sum(ObsProb2(4:6));
            PredProb = [PredProb1, PredProb2Pool];
            ObsProb = [ObsProb1, ObsProb2Pool];
        else % No errors
            PredProb = PredProb1;
            ObsProb = ObsProb1;
        end
        %sum(PredProb)
        %sum(ObsProb)
        PredBuff(i,:) = [PredProb1, PredProb2];
        Ginci = ObsProb .* log(ObsProb ./(PredProb + eps));
        Ginc = sum(Ginci);
        G2 = G2 + Ni(i) * Ginc;
   end;
   G2 = 2 * G2; % Because it is conventional.
   Penalty =  1e3 * (sum(sum(max(P - Pub, 0).^2)) + sum(max(max(Plb - P, 0).^2)));
   if trace
        Penalty
        max(P - Pub, 0)
        max(Plb - P, 0)
        P
   end;
   G2  = G2 + Penalty;   
   AIC = G2 + 2 * sum(Sel);
   BIC = G2 + sum(Sel) * log(ncond * N); 
   %G2
   %sum(Sel)
   %AIC
   %BIC   
end;
%fprintf('ET = %10.3f \n',etime(clock, t0));
end

function [T, G1, G2, Ig1, Ig2, Stat, Qn] = etawiener(P, h, tmax, sw)
% ======================================================================
% Constant drift two-barrier Wiener with drift and starting
% point variability; analytic over eta 
% Usage:
%    [T, G1, G2, Ig1, Ig2, Stat, Qn] = etawiener(P, h, tmax, {switch});
% Parameters:
%    1.  a: upper criterion
%    2.  v: mean drift rate
%    3.  z: starting point
%    4.  eta: drift standard deviation.
%    5.  sz:  starting point standard deviation.
%    6.  Ter: base time.
%    7.  st: Ter variability
% ====================================================================
name = 'ETAWIENER: ';
% format long

if nargin < 4
   sw = 0;
   trace = 0;
end;
if sw == 2
   trace = 1;
   t0 = clock;
else
   trace = 0;
end;
if length(P) ~= 7
   disp([name, ' Incorrect length parameter vector, exiting...']);
   length(P)
   return;
end;
      
epsilon = .0001;
Stat = zeros(1,4);


% Rectangular mass for starting point variability.
n_sz_step =  11;
U = ones(n_sz_step, 1); 
Rmass = U / n_sz_step ; 
Rstep = [-(n_sz_step-1)/2:(n_sz_step-1)/2]' / (n_sz_step-1); 

a = P(1);
mu = P(2);  % v is mu here.
z = P(3);
eta = P(4);
sz = P(5);
Ter = P(6);
st = P(7);
sigma = 0.1; % Ratcliff scale.

if (sz <= epsilon) 
   % -----------------------------
   % No starting point variability
   % -----------------------------
    a1 = a - z;      % a1 := a - z
    a2 = -z;            % a2 := -z            
else 
     %disp('actual')
     %[a, z]
     %disp('transformed')
     Z = z + Rstep * sz;  % <- Rectangular 
     A1 = (a * U - Z);
     A2 = -Z;
     if any(A1 < epsilon) | any(A2 > -epsilon)
        disp('Criterion range error');
        A1'
        A2'
        Rstep * sz
        return
     end
     %if (Z(1) < epsilon | Z(n_sz_step) > a - epsilon)
     %   disp([Z(1),Z(n_sz_step),a]);
     %   disp('Criterion range error');
     %end;                          
end

if (sz <= epsilon)
   % ------------------------------
   % No starting point variability.
   % ------------------------------
   %Pi= [mu, a1, a2, eta, sigma]
   [T, G1, G2, Ig1, Ig2, Stati] = ewiener(mu, a1, a2, eta, sigma, h, tmax);
   Stat(1:2) = Stati(1:2); 
   Stat(3) = Stati(1) * Stati(3);
   Stat(4) = Stati(2) * Stati(4);
else
   % ------------------------------------------
   % Numerical integration over starting point.
   % -------------------------------------------
    G1 = zeros(1,301);
    G2 = zeros(1,301);    
    Ig1 = zeros(1,301);    
    Ig2 = zeros(1,301); 
    for i = 1:n_sz_step
         [T, Gi1, Gi2, Igi1, Igi2, Stati] = ewiener(mu, A1(i), A2(i), eta, sigma, h, tmax);
         G1 = G1 + Rmass(i) * Gi1;
         G2 = G2 + Rmass(i) * Gi2;
         Ig1 = Ig1 + Rmass(i) * Igi1;
         Ig2 = Ig2 + Rmass(i) * Igi2;
         %Stati
         %[mu, A1(i), A2(i)]
         %plot(T, Gi1, T, Gi2)
         %pause
         Stat(1) = Stat(1) + Rmass(i) * Stati(1);
         Stat(2) = Stat(2) + Rmass(i) * Stati(2);
         Stat(3) = Stat(3) + Rmass(i) * Stati(1) * Stati(3);
         Stat(4) = Stat(4) + Rmass(i) * Stati(2) * Stati(4);
    end
    %plot(T, G1/Stat(1), T, G2/Stat(2))
    %pause
end

Stat(3) = Stat(3) / max(Stat(1), epsilon) + Ter;    % E[T1]
Stat(4) = Stat(4) / max(Stat(2), epsilon) + Ter;    % E[T2]
%Stat
%pause
T = T + Ter;
% --------------------
% Convolve with Ter.
% --------------------
if st > 2 * h
   m = round(st/h);
   n = length(T);
   fe = ones(1, m) / m;
   G1 = conv(G1, fe); 
   G2 = conv(G2, fe);
   G1 = G1(1:n);
   G2 = G2(1:n);
   Ig1 = conv(Ig1, fe); 
   Ig2 = conv(Ig2, fe);
   Ig1 = Ig1(1:n);
   Ig2 = Ig2(1:n);
end;

m1 = max(Ig1);
m2 = max(Ig2);
Ig1 = Ig1 / (m1 + m2);
Ig2 = Ig2 / (m1 + m2);
% Calculate quantiles.
if sw
   Qf5 = [.1, .3, .5, .7, .9]; % Summary quantiles (can be changed).

   % Normalise distributions, select subrange spanning quantiles.
   F1 = Ig1 / max(Ig1); 
   F2 = Ig2 / max(Ig2);
   Dummy5 = [.001,.002,.003,.004,.005];

   I1 = (F1 >= .025 & F1 <= .975);  
   if min(diff(F1(I1))) <= 0 
        Q1 = Dummy5;
        disp('Cannot compute F1 quantiles.');    
   else
        Q1=interp1(F1(I1)', T(I1)', Qf5);
   end;

   I2 = (F2 >= .025 & F2 <= .975);  
   if min(diff(F2(I2))) <= 0 
        Q2 = Dummy5;
        disp('Cannot compute F2 quantiles.');    
   else
        Q2=interp1(F2(I2)', T(I2)', Qf5);
   end
   Qn = [Q1',Q2'];
else
   Qn = [];
end
end

function [T, G1, G2, Ig1, Ig2, Stat] = ewiener(mu, a1, a2, eta, sigma, h, tmax)
% ============================================================================
%  First-passage time predictions for a Wiener diffusion proces, with
%  analytic drift variability. (Port of Julia code, August 29, 2017).
%  Usage:
%        [T, G1, G2, Ig1, Ig2, Stat] = ewiener(mu, a1, a2, eta, sigma, h, tmax)
%  Matlab translation of the auto Pascal-to-C code.
% ============================================================================
T = 0:h:tmax;
ssz = length(T);
G1 = zeros(1, ssz);
G2 = zeros(1, ssz);
Ig1 = zeros(1, ssz);
Ig2 = zeros(1, ssz);
Stat = zeros(1,4);
a2 = abs(a2);   
minus_a2 = -a2;

if eta < 0.001
   % Trap underflow
   eta = 0.001;
end
eta_on_sigma2 = (eta/sigma)^2;
mu_on_eta2 = -(mu/eta)^2;
v0a = (a1 * eta_on_sigma2 + mu) / eta;
v0b = (minus_a2 * eta_on_sigma2 + mu) / eta;
w1d = sigma / (a1 + a2);
w3 = pi * w1d^2;
w1e = pi * sigma / (a1 + a2);
w4 = -0.5 * w1e^2;
w5a = pi * a1 / (a1 + a2);
w5b = pi * a2 / (a1 + a2);

ptot = 0.0;
pa = 0.0;
pb = 0.0;
na = 0.0;
nb = 0.0;
G1(1) = 0.0;
G2(1) = 0.0;
Ig1(1) = 0.0;
Ig2(1) = 0.0;
j = 2;
pat = 1.0; % Initialize arbitrarily (was DO loop)
pbt = 1.0;
while (pat + pbt >= 1e-8 || ptot <= 0.9) && j <= ssz
    t = (j - 1) * h;
    pat = 0.0;
    pbt = 0.0;
    k = 1;
    c3 = 1; % Initialize arbitrarily (was DO loop)
    while c3 >= 1e-20
      v1a = mu_on_eta2 + v0a^2 / (t * eta_on_sigma2 + 1.0);
      v1b = mu_on_eta2 + v0b^2 / (t * eta_on_sigma2 + 1.0);
      v2a = w3 * exp(0.5 * v1a);
      v2b = w3 * exp(0.5 * v1b);
      v3a = v2a / sqrt(t * eta_on_sigma2 + 1.0);
      v3b = v2b / sqrt(t * eta_on_sigma2 + 1.0);
      c3 = k * exp(w4 * k^2 * t);
      c4a = sin(k * w5a);
      c4b = sin(k * w5b);
      
      patk = v3a * c3 * c4a;
      pbtk = v3b * c3 * c4b;     

      pat = pat + patk;
      pbt = pbt + pbtk;
      k = k + 1;
    end;
    % Trap small numerical underflows
    if pat > 0 
         G1(j) = pat;
    else 
         G1(j) = 0.0;
    end
    if pbt > 0
         G2(j) = pbt;
    else 
         G2(j) = 0.0;
    end         
    Ig1(j) = Ig1(j-1) + (G1(j) + G1(j-1)) * h / 2.0;
    Ig2(j) = Ig2(j-1) + (G2(j) + G2(j-1)) * h / 2.0;
    pa = pa + pat;
    pb = pb + pbt;
    ptot = ptot + pat + pbt;
    na = na + t * pat;
    nb = nb + t * pbt;
    j = j + 1;
end
 
%  -----------------------------------------------------------------------------
%   If all mass accounted for before ssz, fill arrays with appropriate values.
%  -----------------------------------------------------------------------------
j = j - 1; 
if j < ssz
    for i = j:ssz 
         G1(i) = 0;
         G2(i) = 0;
         Ig1(i) = Ig1(j);
         Ig2(i) = Ig2(j);
    end
end
pa = pa * h;
pb = pb * h;
na = na * h;
nb = nb * h;
Stat(1) = pa; 
Stat(2) = pb; 
if pa > 0
    Stat(3) = na / pa;
else
    Stat(3) = 0.0;
end
if pb > 0
    Stat(4) = nb / pb;
else
    Stat(4) = 0.0;
end  
end




