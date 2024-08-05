%% Compute the fitness function
function c = fitfun_CB(x,bet,eta,k,son,t,dt,lan1,lan2,IdealR3f,IdealV3f)
x(1)=6;
t1=x(1);deltat=x(2);nu=x(3);
% t2 end time of vaccination and t3 is the end time of disease.
t2=t1+deltat;t3=son;
%% Solve SIR until t1
tt1=(0:dt:t1);R1i=0;I1i=1/10^k;S1i=1-I1i-R1i;V1i=0;
fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[0 t1],[S1i I1i R1i V1i]);
T=sol.x';U=sol.y'; UU1=deval(tt1,sol); UU1=UU1';
S1=UU1(:,1);I1=UU1(:,2);R1=UU1(:,3);V1=UU1(:,4);S1f=S1(end,1);I1f=I1(end,1);R1f=R1(end,1);
%% Solve SIR t1 until t2
tt2=(t1:dt:t2);R2i=R1f;I2i=I1f;V2i=0;S2i=S1f;
fun=@(t,y) sira(t,y,bet,eta,nu);sol=ode45(fun,[t1 t2],[S2i I2i R2i V2i]);
T=sol.x';U=sol.y'; UU2=deval(tt2,sol); UU2=UU2';
S2=UU2(:,1);I2=UU2(:,2);R2=UU2(:,3);V2=UU2(:,4);S2f=S2(end,1);I2f=I2(end,1);R2f=R2(end,1);V2f=V2(end,1);
%% Solve SIR t2 until t3
tt3=(t2:dt:t3);R3i=R2f;I3i=I2f;V3i=V2f;S3i=1-R3i-I3i-V3i;
fun=@(t,y) sira(t,y,bet,eta,0);sol=ode45(fun,[t2 t3],[S3i I3i R3i V3i]);
T=sol.x';U=sol.y'; UU3=deval(tt3,sol); UU3=UU3';
S3=UU3(:,1);I3=UU3(:,2);R3=UU3(:,3);V3=UU3(:,4);S3f=S3(end,1);I3f=I3(end,1);R3f=R3(end,1);V3f=V3(end,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=max(lan1*(R3f-IdealR3f), lan2*(V3f-IdealV3f))+0.01*(R3f-IdealR3f+V3f-IdealV3f);
end 