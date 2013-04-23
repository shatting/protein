function [ prediction ] = fromampl(probname,s,newclasses,opt)
%[ prediction ] = fromampl(probname)
%opt = 0: usual fromampl just gets c, t, tp
%opt = 1: extensive fromampl gets all bond vector variables (except bond
%itself)
% s is the secondary structure sequence. needed for geometryz2geometry
if nargin < 4
    opt = 0;
end


global amploutdir;

c = load('-ascii',[amploutdir,probname,'_c.out']);
z = load('-ascii',[amploutdir,probname,'_z.out']);
%tp = load('-ascii',[amploutdir,probname,'_tp.out']);

if opt == 1
    bond1 = load('-ascii',[amploutdir,probname,'_bond1.out']);
    bond2 = load('-ascii',[amploutdir,probname,'_bond2.out']);
    bond3 = load('-ascii',[amploutdir,probname,'_bond3.out']);
    calphas1 = load('-ascii',[amploutdir,probname,'_calphas1.out']);
    calphas2 = load('-ascii',[amploutdir,probname,'_calphas2.out']);
    calphas3 = load('-ascii',[amploutdir,probname,'_calphas3.out']);
  %  bond = load('-ascii',[amploutdir,probname,'_bond.out']);
  %  calphas = load('-ascii',[amploutdir,probname,'_calphas.out']);
    sc = load('-ascii',[amploutdir,probname,'_s.out']);
    phi = load('-ascii',[amploutdir,probname,'_phi.out']);
    alph = load('-ascii',[amploutdir,probname,'_alph.out']);
    ti = load('-ascii',[amploutdir,probname,'_ti.out']);
    tprime = load('-ascii',[amploutdir,probname,'_tprime.out']);
    Vhydro = load('-ascii',[amploutdir,probname,'_Vhydro.out']);
    %bu1 = load('-ascii',[amploutdir,probname,'_bu1.out']);
    %bu2 = load('-ascii',[amploutdir,probname,'_bu2.out']);
    %bu3 = load('-ascii',[amploutdir,probname,'_bu3.out']);
    %a1 = load('-ascii',[amploutdir,probname,'_a1.out']);
    %a2 = load('-ascii',[amploutdir,probname,'_a2.out']);
    %a3 = load('-ascii',[amploutdir,probname,'_a3.out']);
    %deta = load('-ascii',[amploutdir,probname,'_deta.out']);
    %ainv11 = load('-ascii',[amploutdir,probname,'_ainv11.out']);
    %ainv21 = load('-ascii',[amploutdir,probname,'_ainv21.out']);
    %ainv31 = load('-ascii',[amploutdir,probname,'_ainv31.out']);
    %ainv12 = load('-ascii',[amploutdir,probname,'_ainv12.out']);
    %ainv22 = load('-ascii',[amploutdir,probname,'_ainv22.out']);
    %ainv32 = load('-ascii',[amploutdir,probname,'_ainv32.out']);
    %ainv13 = load('-ascii',[amploutdir,probname,'_ainv13.out']);
    %ainv23 = load('-ascii',[amploutdir,probname,'_ainv23.out']);
    %ainv33 = load('-ascii',[amploutdir,probname,'_ainv33.out']);
    %bn1 = load('-ascii',[amploutdir,probname,'_bn1.out']);
    %bn2 = load('-ascii',[amploutdir,probname,'_bn2.out']);
    %bn3 = load('-ascii',[amploutdir,probname,'_bn3.out']);
end %if opt == 2

numc = length(c);

w = c(1:numc);
if newclasses == 0
    prediction.co = w/100;
else 
    prediction.co = w;
end
%prediction.cp = c(2:numc);
prediction.zo = z;
%prediction.tprime = tp;

global alphas;
if newclasses == 0
    prediction.alphaso = [alphas(:,1:2),alphas(:,3:5)*pi/100];
end
prediction.s = s;

prediction.Readme = ['data found by ampl. sc is the s used in geometry2bond, s is the sec. str. class seq.'];

if opt == 1
    prediction.sc = sc;
    %prediction.bu=[bu1,bu2,bu3];
    %prediction.a = [a1,a2,a3];
    %prediction.deta = deta;
    %only last ainv
    %prediction.ainv = [[ainv11(end),ainv12(end),ainv13(end)];[ainv21(end),ainv22(end),ainv23(end)];[ainv31(end),ainv32(end),ainv33(end)]];
    prediction.bond = [bond1,bond2,bond3];
    %prediction.bn = [bn1,bn2,bn3];
    prediction.calphas = [calphas1,calphas2,calphas3];
    prediction.Vhydro = Vhydro;
    %prediction.bond = bond;
    %prediction.calphas = calphas;
    prediction.phi = phi;
    prediction.ti = ti;
    prediction.tprime = tprime;
    prediction.alph = alph;
    %prediction
end %if opt == 2