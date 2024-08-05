%% This is the main code file, in which built in ga() function is used and also
% Objective function  is coded separately and called in the main file.
clc; close all; clear all;input_parameters;
sonsonuc=[];son1=[];son2=[];son3=[];% memory cell to control and save the results in a matrix 
%% Here lan1 is normalized patient care cost and lan2=l-lan1 is vaccination
% normalized cost. The parameters bet(transmission rate), eta(recovery rate) 
% and nu(vaccination) are constants. k is also a constant and will be
% used in fitfun() as I1i=1/10^k. son is ending time of disease. t is time span 
% and dt are step sizes.
for lan1=0.01:0.01:0.99
    lan2=1-lan1;eta=3;bet=6;k=15;son=40;t=(0:1:son)';dt=0.1;
    
% Calling the fitfun() function, as it is our objective function.   
ObjectiveFunction = @(x)fitfun(x,bet,eta,k,son,t,dt,lan1,lan2);
% Number of variables nvar in this case 3. LB and UB are lower and upper bound
% respectively.
nvars = 3;  
% Creation of genetic algorithm options structure to handle all the
% options.
options = gaoptimset('PopulationType', 'doubleVector', ...
    'PopulationSize', 100, ... 
    'EliteCount', '0.05*PopulationSize', ...  
    'CrossoverFraction', 0.8, ...
    'Generations', 300, ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-3, ...
    'CreationFcn',@gacreationlinearfeasible, ...
    'FitnessScalingFcn', @fitscalingrank, ...
    'SelectionFcn', @selectionstochunif, ... 
    'CrossoverFcn',@crossoverarithmetic, ...
    'MutationFcn',@mutationadaptfeasible, ...
    'Display', 'final', ...
    'PlotFcns', @gaplotbestf);

% Calling the ga() function, to get optimum value of objective function.
% finds a local minimum, x, to the objective function also returns fval, 
% the value of the objective function at x. it also returns exitflag, an 
% integer identifying the reason the algorithm terminated, and output, a 
% structure that contains output from each generation and other information 
% about the performance of the algorithm. It also returns a matrix population, 
% whose rows are the final population, and a vector scores, the scores of 
% the final population.
[x,fval,exitFlag,output,population,scores] = ga(ObjectiveFunction,nvars,...
[],[],[],[],LB,UB,[],[],options);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t1 is start of vaccination period deltat is duration of vaccination 
% period and nu is percentage of vaccinated population. These values are
% generated by ga() randomly in an x vector i.e. x=[x(1) x(2) x(3)].
t1=x(1,1);deltat=x(1,2);nnu=x(1,3);t2=t1+deltat;t3=son;sonuc=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve SIR until t1
% setting Initial values and time span (on which integration is required) 
% for ode45() solver.[S1i I1i R1i V1i] are the initial values for
% susceptibles, Infected, Removed(recovered) and vaccinated. Solve the 
% problem with no vaccination 
tt1=(0:dt:t1);
R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;
% Calling the sira() function, that includes the compartmental ordinary
% differential equation system representing SIRV model.
fun=@(t,y) sira(t,y,bet,eta,0);% with no vaccination so nu=0 
% Calling the ode45() function, that will solve the compartmental ordinary
% differential equation system representing SIRV model.
sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
T=sol.x';U=sol.y'; % These are column vectors 
% Calling the deval() function, that will evaluate the solution sol of above 
% differential equation problem at the points contained in tt1.
UU1=deval(tt1,sol); UU1=UU1';
S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);
% [S1f I1f R1f V1f] are the final values for susceptibles, Infected, 
% Removed(recovered) and vaccinated.
S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);
son1=[son1;S1f I1f R1f];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve SIR t1 until t2
% setting Initial values and time span (on which integration is required) 
% for ode45() solver.[S2i=S1f I2i=I1f R2i=R1f V2i=0] are the initial 
% values for susceptibles, Infected, Removed(recovered) and vaccinated. 
tt2=(t1:dt:t2);
R2i=R1f;I2i=I1f;V2i=0;S2i=S1f;
if S1f<eta/bet, nu=0;else nu=nnu; end
% Calling the sira() function, that includes the compartmental ordinary
% differential equation system representing SIRV model.
fun=@(t,y) sira(t,y,bet,eta,nu);
sol=ode45(fun,[t1 t2],[S2i I2i R2i V2i]);
T=sol.x';U=sol.y'; % These are column vectors 
% Calling the deval() function, that will evaluate the solution sol of above 
% differential equation problem at the points contained in tt1.
UU2=deval(tt2,sol); UU2=UU2';
S2=UU2(:,1);I2=UU2(:,2);R2=UU2(:,3);V2=UU2(:,4);
% [S2f I2f R2f V2f] are the final values for susceptibles, Infected, 
% Removed(recovered) and vaccinated.
S2f=S2(end,1);I2f=I2(end,1);R2f=R2(end,1);V2f=V2(end,1);
son2=[son2;S2f I2f R2f V2f];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve SIR t2 until t3
% setting Initial values and time span (on which integration is required) 
% for ode45() solver.[S3i=1-R3i-I3i-V3i I3i=I2f R3i=R2f V3i=V2f] are the initial 
% values for susceptibles, Infected, Removed(recovered) and vaccinated. 
tt3=(t2:dt:t3);
R3i=R2f;I3i=I2f;V3i=V2f;S3i=1-R3i-I3i-V3i;
% Calling the sira() function, that includes the compartmental ordinary
% differential equation system representing SIRV model.
fun=@(t,y) sira(t,y,bet,eta,0);
sol=ode45(fun,[t2 t3],[S3i I3i R3i V3i]);
T=sol.x';U=sol.y'; % These are column vectors 
% Calling the deval() function, that will evaluate the solution sol of above 
% differential equation problem at the points contained in tt1.
UU3=deval(tt3,sol); UU3=UU3';
S3=UU3(:,1);I3=UU3(:,2);R3=UU3(:,3);V3=UU3(:,4);
% [S3f I3f R3f V3f] are the final values for susceptibles, Infected, 
% Removed(recovered) and vaccinated.
S3f=S3(end,1);I3f=I3(end,1);R3f=R3(end,1);V3f=V3(end,1);
son3=[son3;S3f I3f R3f V3f];
% sonsonuc is matrix to store all the results to use them for graphing.
sonsonuc=[sonsonuc;t1 deltat nnu fval S3f I3f R3f V3f lan1 1-lan1];
end


%%Plot between R3f and V3f of GA
figure()
plot(sonsonuc(:,7),sonsonuc(:,8),'*')
title('Rf vs Vf with GA-ODE')
xlabel('Rf')
ylabel('Vf')
sonuct1    =sonsonuc(:,1);
sonucdeltat=sonsonuc(:,2);
sonucnnu   =sonsonuc(:,3);
sonucS3f   =sonsonuc(:,5);
sonucI3f   =sonsonuc(:,6);
sonucR3f   =sonsonuc(:,7);
sonucV3f   =sonsonuc(:,8);

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct computation
eta=3;bet=6;k=15;son=40;t=(0:1:son)';dt=0.1;
% No vaccination  full set
tt1=(0:dt:son);t1=son;
R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;
fun=@(t,y) sira(t,y,bet,eta,0);
sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
T=sol.x';U=sol.y'; % These are column vectors 
UU1=deval(tt1,sol); UU1=UU1';
S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);
S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);
figure(),plot(tt1,UU1)
refRf=R1f;
Sm=eta/bet;
Rm=-eta/bet*log(Sm);
Im=1-Rm-Sm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vaccination
sonuc=[];sonuc1=[];sonuc2= [];
for t1=10:0.5:15
    for deltat=1:0.5:10
        for nnu=0.2:0.01:1
        
% % Solve SIR until t1
t2=t1+deltat;t3=son;
tt1=(0:dt:t1);
R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;
fun=@(t,y) sira(t,y,bet,eta,0);
sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
T=sol.x';U=sol.y'; % These are column vectors 
UU1=deval(tt1,sol); UU1=UU1';
S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);
S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);   
sonuc1=[sonuc1; t1 deltat nnu  S1f I1f R1f];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt2=(t1:dt:t2);
%if S1f<eta/bet, nu=0;else nu=nnu;end
nu=nnu;
R2i=R1f;I2i=I1f;V2i=0;S2i=S1f;
fun=@(t,y) sira(t,y,bet,eta,nu);
sol=ode45(fun,[t1 t2],[S2i I2i R2i V2i]);
T=sol.x';U=sol.y'; % These are column vectors 
UU2=deval(tt2,sol); UU2=UU2';
S2=UU2(:,1);I2=UU2(:,2);R2=UU2(:,3);V2=UU2(:,4);
S2f=S2(end,1);I2f=I2(end,1);R2f=R2(end,1);V2f=V2(end,1);
sonuc2=[sonuc2; t1 deltat nnu  S2f I2f R2f V2f];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt3=(t2:dt:t3);
R3i=R2f;I3i=I2f;V3i=V2f;S3i=1-R3i-I3i-V3i;
fun=@(t,y) sira(t,y,bet,eta,0);
sol=ode45(fun,[t2 t3],[S3i I3i R3i V3i]);
T=sol.x';U=sol.y'; % These are column vectors 
UU3=deval(tt3,sol); UU3=UU3';     
S3=UU3(:,1);I3=UU3(:,2);R3=UU3(:,3);V3=UU3(:,4);
S3f=S3(end,1);I3f=I3(end,1);R3f=R3(end,1);V3f=V3(end,1);
%plot(tt1,UU1,tt2,UU2,tt3,UU3),pause
sonuc=[sonuc;t1 deltat nnu S3f I3f R3f V3f ];
end;end;end

%% Plot between R3f and V3f of direct computation
figure()
plot(sonuc(:,6),sonuc(:,7),'.')
title('Rf vs Vf with Direct Computation')
xlabel('Rf')
ylabel('Vf')
sonuct1    =sonuc(:,1);
sonucdeltat=sonuc(:,2);
sonucnnu   =sonuc(:,3);
sonucS3f   =sonuc(:,4);
sonucI3f   =sonuc(:,5);
sonucR3f   =sonuc(:,6);
sonucV3f   =sonuc(:,7);

%% Plot between t1 and R3f of direct computation
%figure(),plot(sonuc(:,1),sonuc(:,6),'.')
%title('Plot between t1 and R3f of direct computation')
%xlabel('t1')
%ylabel('R3f')

%% Plot between deltat and R3f of direct computation
%figure(),plot(sonuc(:,2),sonuc(:,6),'.')
%title('Plot between deltat and R3f of direct computation')
%xlabel('deltat')
%ylabel('R3f')

%% Plot between nuu and R3f of direct computation
%figure(),plot(sonuc(:,3),sonuc(:,6),'.')
%title('Plot between nuu and R3f of direct computation')
%xlabel('nuu')
%ylabel('R3f')

%% Plot between t1 and V3f of direct computation
%figure(),plot(sonuc(:,1),sonuc(:,7),'.')
%title('Plot between t1 and V3f of direct computation')
%xlabel('t1')
%ylabel('R3f')

%% Plot between deltat and V3f of direct computation
%figure(),plot(sonuc(:,2),sonuc(:,7),'.')
%title('Plot between deltat and V3f of direct computation')
%xlabel('deltat')
%ylabel('R3f')

%% Plot between nuu and V3f of direct computation
%figure(),plot(sonuc(:,3),sonuc(:,7),'.')
%title('Plot between nuu and V3f of direct computation')
%xlabel('nuu')
%ylabel('R3f')
%temp=[];
%% in direct computation case alpha is equivalent to lan1
%for alpha=0.1:0.1:0.9
    % cost calculation
     % temp=[temp alpha*sonucR3f+(1-alpha)*sonucV3f];
    %figure()
    %plot(alpha*sonucR3f+(1-alpha)*sonucV3f),title(alpha)
    %title('Plot alpha*sonucR3f+(1-alpha)*sonucV3f of direct computation with lan1 stepsize=0.1 of direct computation')
    %pause
%end


A=[];
for alpha=0.1:0.1:0.9
    % cost calculation
    z=alpha*(sonucR3f)+(1-alpha)*(sonucV3f);
    %plot(z),axis([0 20000 0 1]);
 [a b]=min(z);
 A=[A;a b alpha];
end

ara=A(:,2);
tun=[];
for i=1:size(ara,1)
    tun=[tun; sonuct1(ara(i,1),1) sonucdeltat(ara(i,1),1) sonucnnu(ara(i,1),1) A(i,2) sonucS3f(ara(i,1),1)  sonucI3f(ara(i,1),1) sonucR3f(ara(i,1),1)  sonucV3f(ara(i,1),1)];
end

%% Plot between minimum values of R3f and V3f of direct computation
figure()
plot(tun(:,7),tun(:,8),'r*')
title('Plot between minimum values of R3f and V3f of direct computation')

%% comparison between minimum(R3f and V3f) with all calculated (R3f and V3f) of direct computation
figure()
plot(sonuc(:,6),sonuc(:,7),'.')
title('comparison between minimum(R3f and V3f) with all calculated (R3f and V3f) of direct computation')
xlabel('R3f')
ylabel('V3f')
hold
plot(tun(:,7),tun(:,8),'r*')   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%comparison between plots(R3f vs V3f) of GA and direct computation
figure()
plot(tun(:,7),tun(:,8),'r*')
title('Rf vs Vf plots with GA-ODE and Direct Computation')
xlabel('Rf')
ylabel('Vf')
% Comparison Plot
hold
plot(sonsonuc(:,7),sonsonuc(:,8),'s')
legend('Direct computation','GA-ODE')



