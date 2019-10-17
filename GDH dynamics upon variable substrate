close all;
clear all;
clc;
%% Initial conditions.....
x0 = [0.1 0.1 0.1 0.1 0.1 0.1 0.088 0.1];
%% Time renge......
tspan = 0:0.05:200;eps=0.001;
time_span=size(tspan);
time_last=time_span(1,2);
l=time_last;
%%constant paramiters....
k7=0.480;k9=5.540;k7d=1.130;k10d=2.800;k1=5.7;k2=1.3669;k3=5.5;k4=2.81;k5=0.025;k6=12;
k8=5.77;k10=0.6964;k11=2.03;k12=1.45*k11;k8d=1.4415;k9d=1.0;
% NH4 input from outside......
% wild type E-coli NH4 value (2 mM) ref: Doucette at. all..Nature chemical biology..
% KG input from outside.......
% Wild tipe E-coli KG value (9.81+-0.62) ref: Doucette at. all..Nature chemical biology..
%% Input parameters......
NADPH=0.10;%KG=0.5;
NADP=0.96;
%% no.1
% wild type............
NH4=10;KG=5;s_index=0;
[T1,X1]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG,k6,k5,NH4,k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
if (X1>zeros(size(X1)))%condition for nonnegativity
    for k=1:10
        if(abs(X1(l+1-k,:)-X1(l-k,:))<eps*ones(size(X1(l-k,:))))%condition for steady state
            s_index=s_index+1; % psudo index for steady-state.
        else
            s_index=0;
        end
    end
else
    disp('system is out of bound:');
end
if(s_index==10)
    disp('system is at steady state');
    s_last1(1,:)= X1(l,:);%steady-state value
else
    s_last1(1,:)=0;
end
%f1 = [NH4 KG NH4/KG s_last1 s_last1(1,5)/s_last1(1,2) s_last1(1,8)/NH4]
xlswrite('table_final1.xls',s_last1);
%%
%%no.2
% NH4 two fold up with variable KG ...
NH4=10; %five fold up
KG=1+49.*rand(1,100);  %KG is variable with 100 uniformly distributed random number
s_index=0;
for i=1:length(KG)
    s_index=0;
    [T2,X2]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG(i),k6,k5,NH4,k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
    if (X2>zeros(size(X2)))%condition for nonnegativity
        for k=1:10
            if(abs(X2(l+1-k,:)-X2(l-k,:))<eps*ones(size(X2(l-k,:))))%condition for steady state
                s_index=s_index+1; % psudo index for steady-state.
            else
                s_index=0;
            end
        end
    else
        disp('system is out of bound:');
    end
    if(s_index==10)
        disp('system is at steady state');
        s_last2(i,:)= X2(l,:);%steady-state value at each KG input
    else
        s_last2(i,:)=0;
    end
end
n_k = NH4./KG;j=1;
M2 = [n_k' s_last2];
for i=1:length(KG)
    if(M2(i,1)>1.5)
        M2_s(j,:) = M2(i,:);
        j=j+1;
    end
end
%f2 = [mean(M2_s(:,9)) std(M2_s(:,9)) mean(M2_s(:,3)) std(M2_s(:,3)) mean(M2_s(:,6)) std(M2_s(:,6)) mean(M2_s(:,6)./M2_s(:,3)) std(M2_s(:,6)./M2_s(:,3)) mean(M2_s(:,9)./NH4) std(M2_s(:,9)./NH4)];
xlswrite('table_final2.xls',M2_s);
%%
%%no.3
% NH4 down to 2 and KG is random variable
NH4=5;
KG=1+49.*rand(1,100);  %KG is variable with 100 uniformly distributed random number
s_index=0;
for i=1:length(KG)
    s_index=0;
    [T3,X3]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG(i),k6,k5,NH4,k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
    if (X3>zeros(size(X3)))%condition for nonnegativity
        for k=1:10
            if(abs(X3(l+1-k,:)-X3(l-k,:))<eps*ones(size(X3(l-k,:))))%condition for steady state
                s_index=s_index+1; % psudo index for steady-state.
            else
                s_index=0;
            end
        end
    else
        disp('system is out of bound:');
    end
    if(s_index==10)
        disp('system is at steady state');
        s_last3(i,:)= X3(l,:);%steady-state value at each KG input
    else
        disp('system is NOT at steady state');
        s_last3(i,:)=0;
    end
end
n_k3 = NH4./KG;j=1;
M3 = [n_k3' s_last3];
for i=1:length(KG)
    if(M3(i,1)>1.5)
        M3_s(j,:) = M3(i,:);
        j=j+1;
    end
end
%f3 = [mean(M3_s(:,9)) std(M3_s(:,9)) mean(M3_s(:,3)) std(M3_s(:,3)) mean(M3_s(:,6)) std(M3_s(:,6)) mean(M3_s(:,6)./M3_s(:,3)) std(M3_s(:,6)./M3_s(:,3)) mean(M3_s(:,9)./NH4) std(M3_s(:,9)./NH4)];
xlswrite('table_final3.xls',M3_s);

%%
%%no.4
%Excess KG with variable NH4
NH4=1+9.0.*rand(1,100);
KG=5;
s_index=0;
for i=1:length(NH4)
    s_index=0;
    [T4,X4]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG,k6,k5,NH4(i),k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
    if (X4>zeros(size(X4)))%condition for nonnegativity
        for k=1:10
            if(abs(X4(l+1-k,:)-X4(l-k,:))<eps*ones(size(X4(l-k,:))))%condition for steady state
                s_index=s_index+1; % psudo index for steady-state.
            else
                s_index=0;
            end
        end
    else
        disp('system is out of bound:');
    end
    if(s_index==10)
        disp('system is at steady state');
        s_last4(i,:)= X4(l,:);%steady-state value at each KG input
    else
        disp('system is not at steady state');
        s_last4(i,:)=0;
    end
end
n_k4 = NH4/KG;j=1;
M4 = [n_k4' s_last4];
for i=1:length(NH4)
    if(M4(i,1)>1.5)
        M4_s(j,:) = M4(i,:);
        j=j+1;
    end
end

%f4 = [mean(M4_s(:,9)) std(M4_s(:,9)) mean(M4_s(:,3)) std(M4_s(:,3)) mean(M4_s(:,6)) std(M4_s(:,6)) mean(M4_s(:,6)./M4_s(:,3)) std(M4_s(:,6)./M4_s(:,3)) mean(M4_s(:,9)./NH4') std(M4_s(:,9)./NH4')];
xlswrite('table_final4.xls',M4_s);


%%
%%no.5
%Excess KG with variable NH4
NH4=1+9.0.*rand(1,100);
KG=2.5;
s_index=0;
for i=1:length(NH4)
    s_index=0;
    [T5,X5]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG,k6,k5,NH4(i),k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
    if (X5>zeros(size(X5)))%condition for nonnegativity
        for k=1:10
            if(abs(X5(l+1-k,:)-X5(l-k,:))<eps*ones(size(X5(l-k,:))))%condition for steady state
                s_index=s_index+1; % psudo index for steady-state.
            else
                s_index=0;
            end
        end
    else
        disp('system is out of bound:');
    end
    if(s_index==10)
        disp('system is at steady state');
        s_last5(i,:)= X5(l,:);%steady-state value at each KG input
    else
        disp('system is not at steady state');
        s_last5(i,:)=0;
    end
end
n_k5 = NH4/KG;j=1;
M5 = [n_k5' s_last5];
for i=1:length(NH4)
    if(M5(i,1)>1.5)
        M5_s(j,:) = M5(i,:);
        j=j+1;
    end
end

%f5 = [mean(M5_s(:,9)) std(M5_s(:,9)) mean(M5_s(:,3)) std(M5_s(:,3)) mean(M5_s(:,6)) std(M5_s(:,6)) mean(M5_s(:,6)./M5_s(:,3)) std(M5_s(:,6)./M5_s(:,3)) mean(M5_s(:,9)./NH4') std(M5_s(:,9)./NH4')];
xlswrite('table_final5.xls',M5_s);

%%
%no.6
% NH4 extremely low and variable KG
NH4=0.01;
KG=1+49*rand(1,100);  %KG is variable with 100 uniformly distributed random number
s_index=0;
for i=1:length(KG)
    s_index=0;
    [T6,X6]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG(i),k6,k5,NH4,k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
    if (X6>zeros(size(X6)))%condition for nonnegativity
        for k=1:10
            if(abs(X6(l+1-k,:)-X6(l-k,:))<eps*ones(size(X6(l-k,:))))%condition for steady state
                s_index=s_index+1; % psudo index for steady-state.
            else
                s_index=0;
            end
        end
    else
        disp('system is out of bound:');
    end
    if(s_index==10)
        disp('system is at steady state');
        s_last6(i,:)= X6(l,:);%steady-state value at each KG input
    else
        disp('system is not at steady state');
        s_last6(i,:)=0;
    end
end
%f6 = [mean(s_last6(:,8)) std(s_last6(:,8)) mean(s_last6(:,2)) std(s_last6(:,2)) mean(s_last6(:,5)) std(s_last6(:,5)) mean(s_last6(:,5)./s_last6(:,2)) std(s_last6(:,5)./s_last6(:,2)) mean(s_last6(:,8)./NH4) std(s_last6(:,8)./NH4)];
xlswrite('table_final6.xls',s_last6);
%%
%%no.7
% NH4 down to 2 and KG is random variable
NH4=1+9.0.*rand(1,100);
KG=100;  %KG is variable with 100 uniformly distributed random number
s_index=0;
for i=1:length(NH4)
    s_index=0;
    [T7,X7]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG,k6,k5,NH4(i),k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
    if (X7>zeros(size(X7)))%condition for nonnegativity
        for k=1:10
            if(abs(X7(l+1-k,:)-X7(l-k,:))<eps*ones(size(X7(l-k,:))))%condition for steady state
                s_index=s_index+1; % psudo index for steady-state.
            else
                s_index=0;
            end
        end
    else
        disp('system is out of bound:');
    end
    if(s_index==10)
        disp('system is at steady state');
        s_last7(i,:)= X7(l,:);%steady-state value at each KG input
    else
        disp('system is NOT at steady state');
        s_last7(i,:)=0;
    end
end
%f7 = [mean(s_last7(:,8)) std(s_last7(:,8)) mean(s_last7(:,2)) std(s_last7(:,2)) mean(s_last7(:,5)) std(s_last7(:,5)) mean(s_last7(:,5)./s_last7(:,2)) std(s_last7(:,5)./s_last7(:,2)) mean(s_last7(:,8)./NH4) std(s_last7(:,8)./NH4)];
xlswrite('table_final7.xls',s_last7);

%%
%disp('NH4 sufficiency:')
%f = [f2;f3;f4;f5]
%disp('NH4 and KG constant:')
%f1
%disp('NH4 limitation:')
%f6
%f7
%%

