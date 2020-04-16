options = optimset(@fminsearch);
options = optimset(options, 'Display', 'iter');
options = optimset(options, 'TolFun', 1e-2);
options = optimset(options, 'TolX', 1e-2);
options = optimset(options, 'FunValCheck', 'on');
%options
disp('Command format is: ');
    disp('pest=fminsearch(@fname, P(Sel==1), options, P(Sel==0), Sel, Data);')
