% The main function for analytic solution for FvCB model 
% first drafted by ZHT 
% Last update by ZHT Mar 26, 2025
% Copyright reserved by ZHT @ YNU

clear all
close all
cd D:\Matlab2022\bin\WorkSpace\ResearchPapers\LeafPhotosynModel\FvCB_Ver5\
odata = readtable('driveVars_test.xlsx'); 
Tsmp = odata.Tsmp;
[D] = rh2vpd(odata.Rh,odata.Ta);
Ta = odata.Ta;
PAR = odata.PAR;
Vc25 = odata.Vc25;
Ca = odata.Ca * 24 / 44;
u = odata.u;
LAI = odata.LAI;
rh = odata.Rh / 100; 

gb = 0.2 * sqrt(u/0.05); 
% gm = 0.25; 
gm = 2000; 

for i = 1:length(Tsmp)
    [params.Vcmax, params.Vj, params.cv, params.cj, params.gamma, params.Rd] = fParameteriz(Vc25(i), Ta(i), PAR(i));
    params.g1 = 9; 
    params.g0 = 0.001; 
    dVar.Q = PAR(i); 
    dVar.T = Ta(i);
    dVar.u = u(i);
    dVar.Ca = Ca(i); 
    dVar.hs = rh(i); 
    dVar.D = D(i); 
    [Aj1(i),Aj2(i),Aj3(i)] = fAnalytic_BBerry_J(params.g0, params.g1, gb(i), gm, params.Rd, params.cj, params.Vj, params.gamma, dVar.hs, dVar.D, dVar.Ca);
    [Ac1(i),Ac2(i),Ac3(i)] = fAnalytic_BBerry_Vc(params.g0, params.g1, gb(i), gm, params.Rd, params.cv, params.Vcmax, params.gamma, dVar.hs, dVar.D, dVar.Ca);
    A(i) = min(Ac3(i),Aj3(i)) - params.Rd; 
    Cs = dVar.Ca - A(i)/gb(i);
    gs(i) = params.g0 + params.g1 * A(i) / Cs * dVar.hs; 
    Cc(i) = dVar.Ca - A(i) * (1/gm + 1/gb(i) + 1/gs(i)); 

end

figure
subplot(121)
plot(Tsmp,Aj1);
hold on
plot(Tsmp,Aj2);
plot(Tsmp,Aj3);
hold off
legend('Aj1','Aj2','Aj3');
subplot(122)
plot(Tsmp,Ac1);
hold on
plot(Tsmp,Ac2);
plot(Tsmp,Ac3);
hold off
legend('Ac1','Ac2','Ac3');
figure
plot(Tsmp,Aj3);
hold on
plot(Tsmp,Ac3);
hold off
legend('Aj3','Ac3');
figure
plot(Tsmp,Cc);
figure
rCcCa = Cc./Ca'; 
plot(Tsmp,rCcCa);
figure
plot(Tsmp,gs);