function osfqpf(Data, Pred, yscale1, yscale2)
% ====================================================================
% QPFs for the OSF data. - Plots by response 
% Usage:
%     osfqpf(Data, Pred, {ymin}, {ymax})
% ====================================================================

% Flip by response
Eix = [5,6,3,4,1,2,11,12,9,10,7,8];
Data(:,8:14) = Data(Eix,8:14);  % Only swap the errors to track response
Pred(:,8:14) = Pred(Eix,8:14);  % Only swap the errors to track response
Data
% Use default or user-supplied y-axis scaling
switch nargin
    case 2
        y1 = 200;
        y2 = 1750;
    case 3
        y1 = 200;
        y2 = yscale1;
    case 4
        y1 = yscale1;
        y2 = yscale2;
    otherwise
        disp('Wrong number of input arguments')
        return;
end;  

axehandle = setfig6;
names = {'Hi.F. Sp.'; 'Hi.F. Ac.'; 'Eq.F. Sp.'; 'Eq.F. Ac'; 'Lo.F. Sp.'; 'Lo.F. Ac'};

qpf(axehandle(1), Data(1:2,:), Pred(1:2,:), y1, 1000, names(1));
qpf(axehandle(2), Data(7:8,:), Pred(7:8,:), y1, 1750, names(2));
qpf(axehandle(3), Data(3:4,:), Pred(3:4,:), y1, 1000, names(3));
qpf(axehandle(4), Data(9:10,:), Pred(9:10,:), y1, 1750, names(4));
qpf(axehandle(5), Data(5:6,:), Pred(5:6,:), y1, 1000, names(5));
qpf(axehandle(6), Data(11:12,:), Pred(11:12,:), y1, 1750, names(6));

end

function qpf(axhandle, Data, Pred, y1, y2, name)
% ====================================================================
% Plots data for one condition of OSF19 Data. Two discriminability
% levels per plot. Handles incomplete (median, missing) data.
% Data and Pred are 2 x 14 
%
% Usage:
%     qpf(Data, Pred, {ymin}, {ymax})
% ymin, ymax are the minimum and maximum scaling on the y-axis.
% If omitted, uses defaults.
% ====================================================================
% ----------------------------
% Millisecond scaling
% ----------------------------
Ix = [3:7,10:14];
Data(:,Ix) = 1000.0*Data(:,Ix);
Pred(:,Ix) = 1000.0*Pred(:,Ix);

if (size(Data) ~= [2,14]) 
   disp('QPF: Data must be a 2 x 14 matrix, exiting...');
   return;
end;

if (size(Pred) ~= [2,14]) 
   disp('QPF: Pred must be a 2 x 14 matrix, exiting...');
   return;
end;      

%axhandle1=setfig1;
sat = 0.6;
red =0.25;
axes(axhandle);
hold on
OProb1 = Data(:,1);
OProb2 = Data(:,8);
OQnt1 =  Data(:,3:7);
OQnt2 =  Data(:,10:14);
PProb = [Pred(:,1);Pred(:,8)];    
PQnt =  [Pred(:,3:7);Pred(:,10:14)];
[~,Px]=sort(PProb);
symbol = ['o', 's', 'd', 'v', '^'];
for i = 1:5
    plot(OProb1, OQnt1(:,i), symbol(i), ...
        'MarkerSize', 6, 'MarkerEdgeColor', [0,sat,0], 'MarkerFaceColor', [0,sat,0]);
    for j = 1:2 
        if OQnt2(j,i) > 0
            plot(OProb2(j), OQnt2(j,i), symbol(i), ...
                'MarkerSize', 6, 'MarkerEdgeColor', [red,0,sat], ...
                'MarkerFaceColor', [red,0,sat]); 
        end
    end           
    plot(PProb(Px), PQnt(Px, i), 'k.-' , 'MarkerSize', 10);
end;
set(gca, 'XLim', [0,1.0]);
set(gca, 'YLim', [y1,y2]);
xlabel('Response Probability')
ylabel('RT Quantile')
label(gca, .45, .85, char(name));
end

function axhandle = setfig1
% ==========================================================================
% setfig1:
% Script to construct a 2 x 1 figure object and set default properties.
% Returns axis handles in axhandle, figure handle in fhandle.
%===========================================================================
fhandle = figure;
pw = 21;  % Reference figure sizes for computing positions.
pl = 29;
set(0,       'ScreenDepth', 1); 
set(fhandle, 'DefaultAxesBox', 'on', ...
             'DefaultAxesLineWidth', 1.5, ...
             'DefaultAxesFontSize', 14, ...
             'DefaultAxesXLim', [0,Inf], ...
             'DefaultAxesYLim', [0,1.0], ...
             'PaperUnits', 'centi', ...
             'PaperType', 'a4', ...
             'PaperPosition', [1, 1, 19, 27], ...
             'Position', [120, 10, 360, 510]);
%  Add these to list ablove to fix axes.
%             'DefaultAxesXLim', [-50,10]
%             'DefaultAxesYLim', [.175,.6]
set(fhandle, 'DefaultLineLineWidth', 1.5, ...
             'DefaultLineColor', [1,1,1], ...
             'DefaultLineLineStyle', '-', ...
             'DefaultLineMarkerSize', 2);
set(fhandle, 'DefaultTextFontSize', 14);
figure(fhandle);
positions =[ 5 15 10 10];
positions(:,1) = positions(:,1) / pw;
positions(:,2) = positions(:,2) / pl;
positions(:,3) = positions(:,3) / pw;
positions(:,4) = positions(:,4) / pl;  % Normalized Units
axhandle=[];
for i=1:1
    axh=axes('Position', positions(i,:));
    axhandle=[axhandle,axh];
end;
end

