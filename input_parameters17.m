% Input parameters
% 2009 first week of Septemper 37. hafta bizde t=0  % incubation+infection is 4.5 days
%eta=1;bet=1.7;k=5.5;son=50;dt=0.1;
eta=7/4.5;R0=1.7;bet=R0*eta;k=5.5;son=50;dt=0.1;
%t=(0:1:son)';
%tt1=(0:dt:son);t1=son;
%R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;
%fun=@(t,y) sira(t,y,bet,eta,0);
%sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
%T=sol.x';U=sol.y'; UU1=deval(tt1,sol); UU1=UU1';
%S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);
%S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);
%figure(),plot(tt1,UU1)
%Sm=eta/bet; Rm=-eta/bet*log(Sm);    Im=1-Rm-Sm;

LB = [5 0.5 0.1]; UB=[15 30  2];