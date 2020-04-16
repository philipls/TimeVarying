function [G2, AIC, BIC, Pred] = osfvwn1x(Pvar, Pfix, Sel, Dmatrix, trace)
% ========================================================================
% Fit m-stage gamma drift model to the OSF19 data using G2. 
% With sigma2, beta free
% Data matrix "Dmatrix" must be 12 x 14 matrix. Uses median error RT
% if quantiles can't be computed.
% Usage:
%    [G2, AIC, BIC, Pred] = osfvw1x(Pvar,Pfix, Sel, Dmatrix, {trace});
%    
%  P = [as, aa, vh, ve, c,  pz,  eta, szs, sza, Ter  st, sigma2, beta, m]
%       1   2   3   4   5   6    7     8   9    10   11   12    13    14
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
name = 'OSFVWN1X: ';
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
np = 14; 
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
sigma2 = P(12);
beta = P(13);
m = P(14);

sigma1 = 0.1;
Pi = [vh+c, beta, m, as, (1+pz)*as/2, eta, szs, Ter, st, sigma1, sigma2; 
      ve+c, beta, m, as, (1+pz)*as/2, eta, szs, Ter, st, sigma1, sigma2;
      vh,   beta, m, as,        as/2, eta, szs, Ter, st, sigma1, sigma2;
      ve,   beta, m, as,        as/2, eta, szs, Ter, st, sigma1, sigma2;
      vh-c, beta, m, as, (1-pz)*as/2, eta, szs, Ter, st, sigma1, sigma2;
      ve-c, beta, m, as, (1-pz)*as/2, eta, szs, Ter, st, sigma1, sigma2;
      vh+c, beta, m, aa, (1+pz)*aa/2, eta, sza, Ter, st, sigma1, sigma2; 
      ve+c, beta, m, aa, (1+pz)*aa/2, eta, sza, Ter, st, sigma1, sigma2;
      vh,   beta, m, aa,        aa/2, eta, sza, Ter, st, sigma1, sigma2;
      ve,   beta, m, aa,      aa/2, eta, sza, Ter, st, sigma1, sigma2;
      vh-c, beta, m, aa, (1-pz)*aa/2, eta, sza, Ter, st, sigma1, sigma2;
      ve-c, beta, m, aa, (1-pz)*aa/2, eta, sza, Ter, st, sigma1, sigma2];
q1min = min(ValidQ1);


ZBounds = [(1-pz)*as/2, (1+pz)*as/2];
ZBounda = [(1-pz)*aa/2, (1+pz)*aa/2];

% Minimum distance of any starting value to either boundary.
szmaxs = 2 * min(ZBounds) - 0.01;
szmaxa = 2 * min(ZBounda) - 0.01;
mina = min(as, aa);

% -----------------------------------------------------------------------------------------------------
%      as   aa    vh       ve     c      pz   eta    szs      sza       Ter     st  sigma2   beta   m -------------------------------------------------------------------------------------------------------
Ub =  [.28, .28,  0.50,  0.65,  0.25    1.0    0.5,  szmaxs,  szmaxa,  q1min,  0.5   0.15      50,   10];
Lb =  [.05, .05, -0.50  -0.65, -0.25   -1.0,     0,      0,        0       0,    0   0.001     5,    1];
Pub = [.25, .25,  0.45,  0.50,  0.20   0.95,  0.45, szmaxs-epsn,szmaxs-epsn, q1min-0.05, .45, 0.1 45, 9.5];
Plb = [.06, .06, -0.45, -0.50, -0.20  -0.95,     0,      0,        0    epsn,   0, 0.01   7, 2];


if any(P - Ub > 0) || any(Lb - P > 0)
   G2 = 1e5 + ...
         1e3 * (sum(max(P - Pub, 0).^2) + sum(max(Plb - P).^2));
   AIC = 1e8;
   BIC = 1e8;
   
   if trace
       P
       Ub
        max(P - Ub, 0)
        max(Lb - P, 0)
   end;         
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
        %Pi(i,:)
        [t,g1,g2,ig1,ig2,StPr,QnPr] = vwien400x(Pi(i,:), tmax, 1);
        %i
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
        PredProb = max(PredProb, 0); % Trap complex values of logarithm
        Ginci = ObsProb .* log(ObsProb ./(PredProb + eps));
        Ginc = sum(Ginci);
        %Ginci
        G2 = G2 + Ni(i) * Ginc;
   end;
   G2 = 2 * G2; % Because it is conventional.
   Penalty =  1e3 * (sum(sum(max(P - Pub, 0).^2)) + sum(max(max(Plb - P, 0).^2)) );
   if trace
        max(P - Pub, 0)
        max(Plb - P, 0)
        P
   end;
   G2  = G2 + Penalty;   
   AIC = G2 + 2 * sum(Sel);
   BIC = G2 + sum(Sel) * log(N);   
end;
%fprintf('ET = %10.3f \n',etime(clock, t0));
end

