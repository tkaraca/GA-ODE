%% This is the main code file, in which built in ga() function is used and also
% Objective function  is coded separately and called in the main file.
clc; close all; clear all;
sonsonuc=[];son1=[];son2=[];son3=[];% memory cell to control and save the results in a matrix 
%% Here lan1 is normalized patient care cost and lan2=l-lan1 is vaccination
% normalized cost. The parameters bet(transmission rate), eta(recovery rate) 
% and nu(vaccination) are constants. k is also a constant and will be
% used in fitfun() as I1i=1/10^k. son is ending time of disease. t is time span 
% and dt are step sizes.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0.1li t1 ve 0.01delta t ve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct computation
eta=1;bet=2;k=8;son=40;t=(0:1:son)';dt=0.1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0.1li t1 ve 0.01delta t ve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct computation
eta=1;bet=1.8;k=6;son=40;t=(0:1:son)';dt=0.1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0.1li t1 ve 0.01delta t ve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct computation
eta=1;bet=1.8;k=5.5;son=32;t=(0:1:son)';dt=0.1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0.1li t1 ve 0.01delta t ve





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%0.01nnudirect cal?smad?
% Vaccination
sonuc=[];sonuc1=[];sonuc2= [];
for t1=5:0.5:15
    for deltat=0.5:0.5:15
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
for alpha=0.01:0.01:0.99
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


