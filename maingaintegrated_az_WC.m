%% This is the main code file, in which built in ga() function is used and also
% Objective function  is coded separately and called in the main file.
clc; close all; clear all;
input_parameters20;nvars = 3;  SONUC=[];
resultstolga=cell(1,10);
%%%%%%benim yaptıklarım
pps=[300] %POPULATİON SİZE
pec=[0.07]%elitecount percentage
pcf=[0.8]  %crossoverfraction
pgn=[150] %generations
x_say=[]
scenario=1
%%%%% tunning için döngü yaptım
for i=1:numel(pps)
    for ii=1:numel(pec)
        for iii=1:numel(pcf)
            for iiii=1:numel(pgn)
      for say=1:10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%yaptıklarım bitti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% orjinal kod
sonsonuc=[];son1=[];son2=[];son3=[];% memory cell to control and save the results in a matrix
% Ideal R3f
lan1=1;    lan2=1-lan1;    t=(0:1:son)';dt=0.1;
% Calling the fitfun() function, as it is our objective function.   
ObjectiveFunction = @(x)fitfun(x,bet,eta,k,son,t,dt,lan1,lan2);
options = gaoptimset('PopulationType', 'doubleVector', ...
    'PopulationSize', pps(i), ... 
    'EliteCount', ceil(pec(ii)*pps(i)), ...  
    'CrossoverFraction', pcf(iii), ...
    'Generations', pgn(iiii), ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-3, ...
    'CreationFcn',@gacreationlinearfeasible, ...
    'FitnessScalingFcn', @fitscalingrank, ...
    'SelectionFcn', @selectionstochunif, ... 
    'CrossoverFcn',@crossoverarithmetic, ...
    'MutationFcn',@mutationadaptfeasible, ...
    'Display', 'final', ...
    'PlotFcns', @gaplotbestf);
[x,fval,exitFlag,output,population,scores] = ga(ObjectiveFunction,nvars,...
[],[],[],[],LB,UB,[],[],options);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%tunning yapıp tüm tabloyu oluşturdum
x(1,1)=6;
x_say=[x_say; x fval pps(i) pec(ii) pcf(iii) pgn(iiii) say scenario] %%%%%%tunnin tablosu

      end
      scenario=scenario+1
            end
        end
    end
end
tunningTableOfIdealR3f=x_say
[a b]=min(x_say(:,4)) %%%minimum fval
b=b(1)
x=(x_say(b,:)) 
bestparametersIdealR3f=x 
bestfvalR3f=x(4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% benim yaptıklarım bitti
t1=x(1,1);deltat=x(1,2);nnu=x(1,3);t2=t1+deltat;t3=son;sonuc=[];
%% Solve SIR until t1
tt1=(0:dt:t1);R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
T=sol.x';U=sol.y'; UU1=deval(tt1,sol); UU1=UU1';
S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);
son1=[son1;S1f I1f R1f];
%% Solve SIR t1 until t2
tt2=(t1:dt:t2);R2i=R1f;I2i=I1f;V2i=0;S2i=S1f;fun=@(t,y) sira(t,y,bet,eta,nnu);sol=ode45(fun,[t1 t2],[S2i I2i R2i V2i]);
T=sol.x';U=sol.y'; UU2=deval(tt2,sol); UU2=UU2';
S2=UU2(:,1);I2=UU2(:,2);R2=UU2(:,3);V2=UU2(:,4);S2f=S2(end,1);I2f=I2(end,1);R2f=R2(end,1);V2f=V2(end,1);
son2=[son2;S2f I2f R2f V2f];
%% Solve SIR t2 until t3
tt3=(t2:dt:t3);R3i=R2f;I3i=I2f;V3i=V2f;S3i=1-R3i-I3i-V3i;fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[t2 t3],[S3i I3i R3i V3i]);
T=sol.x';U=sol.y'; UU3=deval(tt3,sol); UU3=UU3';
S3=UU3(:,1);I3=UU3(:,2);R3=UU3(:,3);V3=UU3(:,4);S3f=S3(end,1);I3f=I3(end,1);R3f=R3(end,1);V3f=V3(end,1);
son3=[son3;S3f I3f R3f V3f];
% sonsonuc is matrix to store all the results to use them for graphing.
%sonsonuc=[sonsonuc;t1 deltat nnu bestfvalR3f S3f I3f R3f V3f lan1 1-lan1];
sonsonucR3f=[t1 deltat nnu bestfvalR3f S3f I3f R3f V3f lan1 1-lan1]
IdealR3f=R3f;

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%benim yaptıklarım
pps=[300] %POPULATİON SİZE
pec=[0.05]%elitecount percentage
pcf=[0.8]  %crossoverfraction
pgn=[150] %generations
x_say=[]
scenario=1
%%%%% tunning için döngü yaptım
for i=1:numel(pps)
    for ii=1:numel(pec)
        for iii=1:numel(pcf)
            for iiii=1:numel(pgn)
      for say=1:10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%yaptıklarım bitti
sonsonuc=[];son1=[];son2=[];son3=[];% memory cell to control and save the results in a matrix
% Ideal V3f
lan1=0;    lan2=1-lan1;    t=(0:1:son)';dt=0.1;
% Calling the fitfun() function, as it is our objective function.   
ObjectiveFunction = @(x)fitfun(x,bet,eta,k,son,t,dt,lan1,lan2);
options = gaoptimset('PopulationType', 'doubleVector', ...
    'PopulationSize', pps(i), ... 
    'EliteCount', ceil(pec(ii)*pps(i)), ...  
    'CrossoverFraction', pcf(iii), ...
    'Generations', pgn(iiii), ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-3, ...
    'CreationFcn',@gacreationlinearfeasible, ...
    'FitnessScalingFcn', @fitscalingrank, ...
    'SelectionFcn', @selectionstochunif, ... 
    'CrossoverFcn',@crossoverarithmetic, ...
    'MutationFcn',@mutationadaptfeasible, ...
    'Display', 'final', ...
    'PlotFcns', @gaplotbestf);
[x,fval,exitFlag,output,population,scores] = ga(ObjectiveFunction,nvars,...
[],[],[],[],LB,UB,[],[],options);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%tunning yapıp tüm tabloyu oluşturdum
x(1,1)=6;
x_say=[x_say; x fval pps(i) pec(ii) pcf(iii) pgn(iiii) say scenario] %%%%%%tunnin tablosu
      end
      scenario=scenario+1
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%benim yaptııklarım
tunningTableOfIdealV3f=x_say
[a b]=min(x_say(:,4)) %%%minimum fval
b=b(1)
x=(x_say(b,:)) 
bestparametersIdealV3f=x
bestfvalV3f=x(4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% benim yaptıklarım bitti
t1=x(1,1);deltat=x(1,2);nnu=x(1,3);t2=t1+deltat;t3=son;sonuc=[];
%% Solve SIR until t1
tt1=(0:dt:t1);R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
T=sol.x';U=sol.y'; UU1=deval(tt1,sol); UU1=UU1';
S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);
son1=[son1;S1f I1f R1f];
%% Solve SIR t1 until t2
tt2=(t1:dt:t2);R2i=R1f;I2i=I1f;V2i=0;S2i=S1f;fun=@(t,y) sira(t,y,bet,eta,nnu);sol=ode45(fun,[t1 t2],[S2i I2i R2i V2i]);
T=sol.x';U=sol.y'; UU2=deval(tt2,sol); UU2=UU2';
S2=UU2(:,1);I2=UU2(:,2);R2=UU2(:,3);V2=UU2(:,4);S2f=S2(end,1);I2f=I2(end,1);R2f=R2(end,1);V2f=V2(end,1);
son2=[son2;S2f I2f R2f V2f];
%% Solve SIR t2 until t3
tt3=(t2:dt:t3);R3i=R2f;I3i=I2f;V3i=V2f;S3i=1-R3i-I3i-V3i;fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[t2 t3],[S3i I3i R3i V3i]);
T=sol.x';U=sol.y'; UU3=deval(tt3,sol); UU3=UU3';
S3=UU3(:,1);I3=UU3(:,2);R3=UU3(:,3);V3=UU3(:,4);S3f=S3(end,1);I3f=I3(end,1);R3f=R3(end,1);V3f=V3(end,1);
son3=[son3;S3f I3f R3f V3f];
% sonsonuc is matrix to store all the results to use them for graphing.
%sonsonuc=[sonsonuc;t1 deltat nnu bestfvalV3f S3f I3f R3f V3f lan1 1-lan1];
sonsonucV3f=[t1 deltat nnu bestfvalV3f S3f I3f R3f V3f lan1 1-lan1];
IdealV3f=V3f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tunningTableOfobj=cell(1,10)
lanIteration=0%lan1 tunningleri takip etmek için koydum
bestparametersOfobj=[]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sonsonuc=[];son1=[];son2=[];son3=[];% memory cell to control and save the results in a matrix
for lan1=0.1:0.1:0.9,    lan2=1-lan1;    t=(0:1:son)';dt=0.1;
    %%%%%%benim yaptıklarım
pps=[300] %POPULATİON SİZE
pec=[0.03]%elitecount percentage
pcf=[0.8]  %crossoverfraction
pgn=[50] %generations
    lanIteration=lanIteration+1
x_say=[]
scenario=1
%%%%% tunning için döngü yaptım
for i=1:numel(pps)
    for ii=1:numel(pec)
        for iii=1:numel(pcf)
            for iiii=1:numel(pgn)
      for say=1:10
% Calling the fitfun() function, as it is our objective function.   
ObjectiveFunction = @(x)fitfun_CB(x,bet,eta,k,son,t,dt,lan1,lan2,IdealR3f,IdealV3f);
options = gaoptimset('PopulationType', 'doubleVector', ...
    'PopulationSize', pps(i), ... 
    'EliteCount', ceil(pec(ii)*pps(i)), ...  
    'CrossoverFraction', pcf(iii), ...
    'Generations', pgn(iiii), ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-3, ...
    'CreationFcn',@gacreationlinearfeasible, ...
    'FitnessScalingFcn', @fitscalingrank, ...
    'SelectionFcn', @selectionstochunif, ... 
    'CrossoverFcn',@crossoverarithmetic, ...
    'MutationFcn',@mutationadaptfeasible, ...
    'Display', 'final',...
    'PlotFcns', @gaplotbestf);
[x,fval,exitFlag,output,population,scores] = ga(ObjectiveFunction,nvars,...
[],[],[],[],LB,UB,[],[],options);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%tunning yapıp tüm tabloyu oluşturdum
x(1,1)=6;
x_say=[x_say; x fval pps(i) pec(ii) pcf(iii) pgn(iiii) say scenario lan1] %%%%%%tunnin tablosu
      end
      scenario=scenario+1
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%benim yaptııklarım
tunningTableOfobj{lanIteration}=x_say
[a b]=min(x_say(:,4)) %%%minimum fval
b=b(1)
x=(x_say(b,:)) 
bestparametersOfobj=[bestparametersOfobj;x]
bestfval=x(4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% benim yaptıklarım bitti

t1=x(1,1);deltat=x(1,2);nnu=x(1,3);t2=t1+deltat;t3=son;sonuc=[];
%% Solve SIR until t1
tt1=(0:dt:t1);R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
T=sol.x';U=sol.y'; UU1=deval(tt1,sol); UU1=UU1';
S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);
son1=[son1;S1f I1f R1f];
%% Solve SIR t1 until t2
tt2=(t1:dt:t2);R2i=R1f;I2i=I1f;V2i=0;S2i=S1f;fun=@(t,y) sira(t,y,bet,eta,nnu);sol=ode45(fun,[t1 t2],[S2i I2i R2i V2i]);
T=sol.x';U=sol.y'; UU2=deval(tt2,sol); UU2=UU2';
S2=UU2(:,1);I2=UU2(:,2);R2=UU2(:,3);V2=UU2(:,4);S2f=S2(end,1);I2f=I2(end,1);R2f=R2(end,1);V2f=V2(end,1);
son2=[son2;S2f I2f R2f V2f];
%% Solve SIR t2 until t3
tt3=(t2:dt:t3);R3i=R2f;I3i=I2f;V3i=V2f;S3i=1-R3i-I3i-V3i;fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[t2 t3],[S3i I3i R3i V3i]);
T=sol.x';U=sol.y'; UU3=deval(tt3,sol); UU3=UU3';
S3=UU3(:,1);I3=UU3(:,2);R3=UU3(:,3);V3=UU3(:,4);S3f=S3(end,1);I3f=I3(end,1);R3f=R3(end,1);V3f=V3(end,1);
son3=[son3;S3f I3f R3f V3f];
% sonsonuc is matrix to store all the results to use them for graphing.
sonsonuc=[sonsonuc;t1 deltat nnu bestfval S3f I3f R3f V3f lan1 1-lan1]
end
% SONUC(say,:,:)=sonsonuc;
% resultstolga{say}=sonsonuc

save('result_t16_R20.mat','sonsonuc')
filename = 'tümR0_t16_sonuclar_tablo.xlsx'
writematrix(sonsonuc,filename,'Sheet',9)
% 
% %%Plot between R3f and V3f of GA
%   S=squeeze(SONUC(1,:,:));
%  figure(),plot(S(:,7),S(:,8),'.'),title('Rf vs Vf with GA-ODE'),xlabel('Rf'),ylabel('Vf')
% hold
% S=squeeze(SONUC(2,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(3,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(4,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(5,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(6,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(7,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(8,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(9,:,:));plot(S(:,7),S(:,8),'*')
% S=squeeze(SONUC(10,:,:));plot(S(:,7),S(:,8),'*')
% 
% print('FSfigure12','-dpng','-r0')


%sonuct1    =sonsonuc(:,1);
%sonucdeltat=sonsonuc(:,2);
%sonucnnu   =sonsonuc(:,3);
%sonucS3f   =sonsonuc(:,5);
%sonucI3f   =sonsonuc(:,6);
%sonucR3f   =sonsonuc(:,7);
%sonucV3f   =sonsonuc(:,8);

%% Plot between t1 and R3f of GA
%figure(),plot(sonuct1,sonucR3f,'.')
%title('Plot between t1 and R3f of GA')
%xlabel('t1')
%ylabel('R3f')



%% Plot between deltat and R3f of GA
%figure(),plot(sonucdeltat,sonucR3f,'.')
%title('Plot between deltat and R3f of GA')
%xlabel('deltat')
%ylabel('R3f')

%% Plot between nuu and R3f of GA
%figure(),plot(sonucnnu,sonucR3f,'.')
%title('Plot between nuu and R3f of GA')
%xlabel('nuu')
%ylabel('R3f')

%% Plot between t1 and V3f of GA
%figure(),plot(sonuct1,sonucV3f,'.')
%title('Plot between t1 and V3f of GA')
%xlabel('t1')
%ylabel('R3f')

%% Plot between deltat and V3f of GA
%figure(),plot(sonucdeltat,sonucV3f,'.')
%title('Plot between deltat and V3f of GA')
%xlabel('deltat')
%ylabel('R3f')

%% Plot between nuu and V3f of GA
%figure(),plot(sonucnnu,sonucV3f,'.')
%title('Plot between nuu and V3f of GA')
%xlabel('nuu')
%ylabel('R3f')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bunları şimdilik kaldırıyorum
% %%%%R3f tablo haline getir
% obvFunValue=tunningTableOfIdealR3f(:,4);
% PopulationSize=tunningTableOfIdealR3f(:,5);
% elitecountPercentage=tunningTableOfIdealR3f(:,6);
% crossoverFraction=tunningTableOfIdealR3f(:,7);
% generations=tunningTableOfIdealR3f(:,8);
% iterationNumber=tunningTableOfIdealR3f(:,9);
% Scenario=tunningTableOfIdealR3f(:,10);
% tblR3f = table(obvFunValue,PopulationSize,elitecountPercentage,crossoverFraction,generations,iterationNumber,Scenario);
% 
% %%%%%%%%% tüm senaryolar
% tblscenarioR3f=grpstats(tblR3f,"Scenario")
% 
% %%%%% her senarryoya ait istatistikler
% tblstatsR3f= grpstats(tblR3f,"Scenario",["mean","var","min","max","meanci","predci"],"DataVars","obvFunValue","Alpha",0.05)
% 
% %%%%overall istatistikler
% tblstatsR3fWhole= grpstats(tblR3f,[],["mean","var","min","max","meanci","predci"],"DataVars","obvFunValue","Alpha",0.05)
% 
% %senaryo ortalamaları ve CI grafiği alpha=0.05
% grpstats(obvFunValue,Scenario,0.05)
% legend("obvFunValue")
% 
% % figure
% % boxplot(obvFunValue,{PopulationSize,elitecountPercentage,crossoverFraction,generations})
% % title('obvFunValue, Grouped by PopulationSize,elitecountPercentage,crossoverFraction,generations')
% 
% %%%senaryo boxplot grafiği 
% figure
% boxplot(obvFunValue,Scenario)
% title('obvFunValue, Scenario')
% 
% 
% %%%%% N-WAY ANOVA
% %anovan(obvFunValue,{PopulationSize,elitecountPercentage,crossoverFraction,generations},"Continuous",3, ...
%     "Varnames",["PopulationSize","elitecountPercentage","crossoverFraction","generations"]);
% 
% %%%%%% N-WAY ANOVA INTERACTION
% varnames = {'PopulationSize','elitecountPercentage','crossoverFraction','generations'}
% [~,~,stats]=anovan(obvFunValue,{PopulationSize,elitecountPercentage,crossoverFraction,generations},4,3, ...
%    varnames)
% 
% %%%%%% N-WAY ANOVA INTERACTION
% varnames = {'PopulationSize','crossoverFraction',}
% [~,~,stats]=anovan(obvFunValue,{PopulationSize,crossoverFraction},1,2, ...
%    varnames)
% 
% %%%%%% FACTOR COMPARISON
% [results,~,~,gnames] = multcompare(stats,'Dimension',[1,2]);
% 
% tblR3fComp = array2table(results,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% tblR3fComp.("Group A") = gnames(tblR3fComp.("Group A"));
% tblR3fComp.("Group B") = gnames(tblR3fComp.("Group B"))
%  %%%Write the tables and results in excell ParameterAnalysisOfIdealR3f.xlsx
%  %%%for IDEAL R3f
% filename = 'ParameterAnalysisOfIdealR3f.xlsx'
% writematrix(tunningTableOfIdealR3f,filename,'Sheet',1)
% writematrix(bestparametersIdealR3f,filename,'Sheet',2)
% writetable(tblscenarioR3f,filename,'Sheet',3)
% writetable(tblstatsR3fWhole,filename,'Sheet',4)
% writetable(tblstatsR3f,filename,'Sheet',5)
% writetable(tblR3fComp,filename,'Sheet',6)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%V3f tablo haline getir
% obvFunValue=tunningTableOfIdealV3f(:,4);
% PopulationSize=tunningTableOfIdealV3f(:,5);
% elitecountPercentage=tunningTableOfIdealV3f(:,6);
% crossoverFraction=tunningTableOfIdealV3f(:,7);
% generations=tunningTableOfIdealV3f(:,8);
% iterationNumber=tunningTableOfIdealV3f(:,9);
% Scenario=tunningTableOfIdealV3f(:,10);
% tblV3f = table(obvFunValue,PopulationSize,elitecountPercentage,crossoverFraction,generations,iterationNumber,Scenario);
% 
% %%%%%%%%% tüm senaryolar
% tblscenarioV3f=grpstats(tblV3f,"Scenario")
% 
% %%%%% her senarryoya ait istatistikler
% tblstatsV3f= grpstats(tblV3f,"Scenario",["mean","var","min","max","meanci","predci"],"DataVars","obvFunValue","Alpha",0.05)
% 
% %%%overall istatistikler
% tblstatsV3fWhole= grpstats(tblV3f,[],["mean","var","min","max","meanci","predci"],"DataVars","obvFunValue","Alpha",0.05)
% 
% %senaryo ortalamaları ve CI grafiği alpha=0.05
% grpstats(obvFunValue,Scenario,0.05)
% legend("obvFunValue")
% 
% %%%senaryo boxplot grafiği 
% figure 
% boxplot(obvFunValue,Scenario)
% title('obvFunValue, Scenario')
% 
% %%%%% N-WAY ANOVA
% %anovan(obvFunValue,{PopulationSize,elitecountPercentage,crossoverFraction,generations},"Continuous",3, ...
%     "Varnames",["PopulationSize","elitecountPercentage","crossoverFraction","generations"]);
% 
% %%%%%% N-WAY ANOVA INTERACTION
% varnames = {'PopulationSize','elitecountPercentage','crossoverFraction','generations'}
% [~,~,stats]=anovan(obvFunValue,{PopulationSize,elitecountPercentage,crossoverFraction,generations},4,3, ...
%    varnames)
% 
% %%%%%% N-WAY ANOVA INTERACTION
% varnames = {'PopulationSize','crossoverFraction',}
% [~,~,stats]=anovan(obvFunValue,{PopulationSize,crossoverFraction},2,3, ...
%    varnames)
% %%%%%% FACTOR COMPARISON
% [results,~,~,gnames] = multcompare(stats,'Dimension',[1,2]);
% 
% tblV3fcomp = array2table(results,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% tblV3fcomp.("Group A") = gnames(tblV3fcomp.("Group A"));
% tblV3fcomp.("Group B") = gnames(tblV3fcomp.("Group B"))
% 
% 
% %%%Write the tables and results in excell ParameterAnalysisOfIdealR3f.xlsx
%  %%%for IDEAL V3f
% filename = 'ParameterAnalysisOfIdealV3f.xlsx'
% writematrix(tunningTableOfIdealV3f,filename,'Sheet',1)
% writematrix(bestparametersIdealV3f,filename,'Sheet',2)
% writetable(tblscenarioV3f,filename,'Sheet',3)
% writetable(tblstatsV3fWhole,filename,'Sheet',4)
% writetable(tblstatsV3f,filename,'Sheet',5)
% writetable(tblV3fcomp,filename,'Sheet',6)
% 
% 
% %%%%bi-obj tablo haline getir
% tunningTableOfIdealbiobj=tunningTableOfobj{1}
% obvFunValue=tunningTableOfIdealbiobj(:,4);
% PopulationSize=tunningTableOfIdealbiobj(:,5);
% elitecountPercentage=tunningTableOfIdealbiobj(:,6);
% crossoverFraction=tunningTableOfIdealbiobj(:,7);
% generations=tunningTableOfIdealbiobj(:,8);
% iterationNumber=tunningTableOfIdealbiobj(:,9);
% Scenario=tunningTableOfIdealbiobj(:,10);
% tblbiobj = table(obvFunValue,PopulationSize,elitecountPercentage,crossoverFraction,generations,iterationNumber,Scenario);
% 
% %%%%%%%%% tüm senaryolar
% tblscenariobiobj=grpstats(tblbiobj ,"Scenario")
% %%%%% her senarryoya ait istatistikler
% tblstatsbiobj= grpstats(tblbiobj,"Scenario",["mean","var","min","max","meanci","predci"],"DataVars","obvFunValue","Alpha",0.05)
% %%%overall istatistikler
% tblstatsbiobjWhole= grpstats(tblbiobj,[],["mean","var","min","max","meanci","predci"],"DataVars","obvFunValue","Alpha",0.05)
% %%%senaryo ortalamaları ve CI grafiği alpha=0.05
% grpstats(obvFunValue,Scenario,0.05)
% legend("obvFunValue")
% %%%senaryo boxplot grafiği 
% figure
% boxplot(obvFunValue,Scenario)
% title('obvFunValue, Scenario')
% %%%%% N-WAY ANOVA
% %anovan(obvFunValue,{PopulationSize,elitecountPercentage,crossoverFraction,generations},"Continuous",3, ...
%     "Varnames",["PopulationSize","elitecountPercentage","crossoverFraction","generations"]);
% 
% %%%%%% N-WAY ANOVA INTERACTION
% varnames = {'PopulationSize','elitecountPercentage','crossoverFraction','generations'}
% [~,~,stats]=anovan(obvFunValue,{PopulationSize,elitecountPercentage,crossoverFraction,generations},4,3, ...
%    varnames)
% 
% %%%%%% N-WAY ANOVA INTERACTION
% varnames = {'PopulationSize','crossoverFraction',}
% [~,~,stats]=anovan(obvFunValue,{PopulationSize,crossoverFraction},1,2, ...
%    varnames)
% 
% %%%%%% FACTOR COMPARISON
% [results,~,~,gnames] = multcompare(stats,'Dimension',[1,2]);
% 
% tblzComp = array2table(results,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% tblzComp.("Group A") = gnames(tblzComp.("Group A"));
% tblzComp.("Group B") = gnames(tblzComp.("Group B"))
% %%%Write the tables and results in excell ParameterAnalysisOfIdealR3f.xlsx
%  %%%for IDEAL V3f
% filename = 'ParameterAnalysisOfBiObjFun.xlsx'
% writematrix(tunningTableOfIdealbiobj,filename,'Sheet',1)
% writematrix(bestparametersOfobj,filename,'Sheet',2)
% writetable(tblscenariobiobj,filename,'Sheet',3)
% writetable(tblstatsbiobjWhole,filename,'Sheet',4)
% writetable(tblstatsbiobj,filename,'Sheet',5)
% writetable(tblzComp,filename,'Sheet',6)
