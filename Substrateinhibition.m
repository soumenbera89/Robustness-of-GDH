% Substrate inhibition by alpha ketoglutarate......
% Soumen Bera and Amit Chakraborty.......
% E-mail- soumenmath4189@gmail.com.......
% First curve is drown without NADP+ and then 0.96 mM NADP+ is added which
% represent the second curve. NH4 and Alpha KG is inversly varied from 1 to
% 10 mM and 1 to 50 mM consequitively.
clear all;
close all;
clc;
%% Initial conditions.....
x0 = [0.1 .1 .1 .1 .1 .1 .088 .1] ;
%% Time renge......
tspan = 0:0.05:200;
%%constant paramiters....
k7=0.480;k9=5.540;k7d=1.130;k10d=2.800;k1=5.7;k2=1.3669;k3=5.5;k4=2.81;k5=0.025;k6=12;
k8=5.77;k10=0.6964;k11=2.03;k12=1.45*k11;k8d=1.4415;k9d=1.0;
%% Input paramiters......
NADPH=0.10;
NADP=0.0;
eps=0.001;
% Create ammonia gradient.......
NH4_high=10.0;
NH4_low=1.0;
NH4=linspace(NH4_high,NH4_low,100);
% Create a alpha ketoglutarate gradient.....
KG_high=50.0;
KG_low=1.0;
KG=linspace(KG_low,KG_high,100);
% last value calculation from every step iteration.......
time_span=size(tspan);
time_last=time_span(1,2);
l=time_last;
s_last=zeros(length(KG),8);
% Solve differtial equation by RK4 method keeping all solution in positive
% quadrant and establish steady state..
    for i=1:length(KG)
        [T,X]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG(i),k6,k5,NH4(i),k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
        if (X>zeros(size(X)))%condition for nonnegativity
            for k=1:10
                if(abs(X(l+1-k,:)-X(l-k,:))<eps*ones(size(X(l-k,:))))%condition for steady state
                    for j=1:8
                    s_last(i,j)=X(time_last,j);% to take last value for each input value kG
                    end
                else
                    for j=1:8
                    s_last(i,j)=0;
                    end
                end
            end
        else
            disp('system is out of bound:');
        end
    end
for i=1:length(KG)
    const_ratio(i) = s_last(i,5)/s_last(i,2); % calculation of [GDH.NADP+.KG] and [GDH.NADPH.KG] ratio
end
eps1=0.01;
mean_const_ratio = mean(const_ratio);
%...................
% kinetic complex ratio related to alpha-ketoglutarate inhibition
s(1)=1;
for i=1:length(KG)
    if(abs(const_ratio(i)-mean_const_ratio)< eps1)
        s(i)=1;%pusedo index
        disp('result is satisfied \n');
        fprintf('\nconstant ratio=%f \n', const_ratio(i));
        %fprintf('\n input-output ratio\n', s_last(:,8)./KG')
        %fprintf('NH4 = %f \t m_specific = %f\n\n',NH4(i),m_specific(i));
    else
        s(i)=0;
        disp('result is NOT satisfied \n');
    end
end
plot(1./KG, 1./s_last(:,6),'ko--','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);
niceplot;
hold on;
% When NADP+ is added to the mixtute.......
NADP=0.96;
    for i=1:length(KG)
        [T,X]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG(i),k6,k5,NH4(i),k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
        if (X>zeros(size(X)))%condition for nonnegativity
            for k=1:10
                if(abs(X(l+1-k,:)-X(l-k,:))<eps*ones(size(X(l-k,:))))%condition for steady state
                    for j=1:8
                    s_last(i,j)=X(time_last,j);% to take last value for each input value kG
                    end
                else
                    for j=1:8
                    s_last(i,j)=0;
                    end
                end
            end
        else
            disp('system is out of bound:');
        end
    end
for i=1:length(KG)
    const_ratio(i) = s_last(i,5)/s_last(i,2); % calculation of [GDH.NADP+.KG] and [GDH.NADPH.KG] ratio
end
eps1=0.01;
mean_const_ratio = mean(const_ratio);
%...................
% kinetic complex ratio related to alpha-ketoglutarate inhibition
s(1)=1;
for i=1:length(KG)
    if(abs(const_ratio(i)-mean_const_ratio)< eps1)
        s(i)=1;%pusedo index
        disp('result is satisfied \n');
        fprintf('\nconstant ratio=%f \n', const_ratio(i));
        %fprintf('\n input-output ratio\n', s_last(:,8)./KG')
        %fprintf('NH4 = %f \t m_specific = %f\n\n',NH4(i),m_specific(i));
    else
        s(i)=0;
        disp('result is NOT satisfied \n');
    end
end
plot(1./KG, 1./s_last(:,6),'ko--','MarkerSize',6);
niceplot;
%print('kg_inhibition','-depsc','-tiff')
saveas(gca, 'kf_inhibition.eps','epsc');
