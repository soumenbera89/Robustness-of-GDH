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
%NADPH=linspace(0.019,0.191,400);
NADPH=0.10;
NADP=0.0;
eps = 0.001;
% Create ammonia gradient.......
NH4_high=10;
NH4_low=1;
NH4=linspace(NH4_high,NH4_low,100);
%NH4=1.0+9.0.*rand(1,100);
%NH4 = 5.0;
KG_high=50.0;
KG_low=1.0;
KG=linspace(KG_low,KG_high,100);
%KG=1+49.0.*rand(1,100);
%KG =10;
%%NH4 constant input.....
%NH4=2;
%KG=10;
time_span=size(tspan);
time_last=time_span(1,2);
l=time_last;
s_last=zeros(length(KG),8);
%% ode solver..
%for i=1:length(NH4)
%for p=1:20
    %k1=a+(b-a)*rand(1);k2=a+(b-a)*rand(1);k3=a+(b-a)*rand(1);k4=a+(b-a)*rand(1);k5=a+(b-a)*rand(1);k6=a+(b-a)*rand(1);
    %k8=a+(b-a)*rand(1);k10=a+(b-a)*rand(1);k11=a+(b-a)*rand(1);k12=1.45*k11;k8d=a+(b-a)*rand(1);k9d=a+(b-a)*rand(1);
    %[T,X]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG(k),k6,k5,NH4(k),k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
    %if(abs(X(time_last,5)-X(time_last,2))<0.5)
     %   k11(p)=k1;k22(p)=k2;k33(p)=k3;k44(p)=k4;k55(p)=k5;k66(p)=k6;k88(p)=k8;k100(p)=k10;k110(p)=k11;k120(p)=k12;
      %  k8dd(p)=k8d;k9dd(p)=k9d;
    %else
     %   k11(p)=0;k22(p)=0;k33(p)=0;k44(p)=0;k55(p)=0;k66(p)=0;k88(p)=0;k100(p)=0;k110(p)=0;k120(p)=0;
      %  k8dd(p)=0;k9dd(p)=0;
       % print('the experiment is fail\n\n');
%for kp=1:length(k5)
    for i=1:length(KG)
        [T,X]=ode45(@(t,x) GDH_system(t,x,k1,NADPH,k4,k2,k3,KG(i),k6,k5,NH4(i),k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d),tspan,x0);
        if (X>zeros(size(X)))%condition for nonnegativity
            for k=1:10
                if(abs(X(l+1-k,:)-X(l-k,:))<eps*ones(size(X(l-k,:))))%condition for steady state
                    for j=1:8
                    s_last(i,j)=X(time_last,j);% to take last value for each input value kG
                    %Etot{i}=X(:,1)+X(:,2)+X(:,3)+X(:,4)+X(:,5)+X(:,6)+X(:,7);
                    %specific{i}=abs((k7d*X(:,6)-k8*X(:,8).*X(:,4)-k8d*X(:,8).*X(:,7)));
                    %m_specific(i)=mean(specific{i}(1:4000));
                    end
                    %KG_inhibition(i)=abs(k7d*X(time_last,6)-k8*X(time_last,8).*X(time_last,4)-k8d*X(time_last,8).*X(time_last,7));
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
    
    % kinetic complex ratio related to alpha-ketoglutarate inhibition
for i=1:length(KG)
    c_ratio(i) = s_last(i,5)/s_last(i,2);
end
eps1=0.1;
f_d(1)=0;
c_const(1)=0;
for i=2:length(KG)
    f_d(i)=abs((c_ratio(i)-c_ratio(i-1))/(NH4(i)-NH4(i-1)));
    if(f_d(i)<eps1)
        c_const(i)=1;
    else
        c_const(i)=0;
    end
end
%...................


% product-input ratio-a link
eps3 = 0.1;
for i=1:length(KG)
    io_ratio(i) = s_last(i,8)/NH4(i);
end
ff_d(1)=0;
io_const(1)=0;
for i=2:length(KG)
    ff_d(i)=abs((io_ratio(i)-io_ratio(i-1))/(NH4(i)-NH4(i-1)));
    if(ff_d(i)<eps3)
        io_const(i)=1;
    else
        io_const(i)=0;
    end
end

robust_matrix(1,:)= [1 0 0 0 0 0 0 0 0 0 0 0];
for i=1:length(KG)
    if(c_const(i)==1 && io_const(i)== 1)
        robust_index(i)=1;
        disp('\n Result is robust \n');
        robust_matrix(i,:)= [i NADP NADPH NH4(i) KG(i) s_last(i,5) s_last(i,2) s_last(i,6) c_const(i) io_const(i) c_ratio(i) io_ratio(i)];
    else
        robust_index(i)=0;
        robust_matrix(i,:)= [i NADP NADPH NH4(i) KG(i) s_last(i,5) s_last(i,2) s_last(i,6) c_const(i) io_const(i) c_ratio(i) io_ratio(i)];
        disp('Robust result is NOT established \n');
    end
end
robust_matrix

%for k=1:19
 %   m_change(k)=m_specific(k+1)-m_specific(k);
%end
%h=figure;
%subplot(1,2,1);
%plot(NH4,m_specific,'ko--');
%plot(NH4,m_specific,'ko--','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
%plot(NH4,m_specific,'o','MarkerEdgeColor','k','MarkerSize',8);
%hold on;
%plot(1./KG,1./s_last(:,8),'o-');
%hold on;
%figure
%plot(KG, s_last(:,6),'ko--','MarkerFaceColor','r','MarkerSize',8);
%hold on;
%plot(KG,1./s_last(:,8));
%hold on;
%plot(s_last(:,5),s_last(:,6),'o-');
%end
%ylabel('Specific Activity of GDH','FontWeight','Bold','FontSize',10);
%xlabel('Extracellular [NH^{+}_{4}]','FontWeight','Bold','FontSize',10);
%niceplot;
%A=[NH4' KG' s_last(:,8) s_last(:,5)];
%xlswrite('KG_inhibition',A);
%%
%Robustness plot
%n=0;
%mm1(1)=0;mm2(1)=0;
%for i=2:length(KG)
 %   if(robust_index(i)==1)
  %      n=n+1;
   %     mm1(i)=c_ratio(i);
    %    mm2(i)=io_ratio(i);
    %else
     %   mm1(i)=0;
      %  mm2(i)=0;
    %end
%end
%m1=sum(mm1)/n;m2=sum(mm2)/n; %mean value
%vv1=0;vv2=0;
%for i=2:length(KG)
 %   if(mm1(i)~=0 && mm2(i)~=0)
  %      vv1=vv1+(m1-mm1(i))^2;
   %     vv2=vv2+(m2-mm2(i))^2;
    %end
%end
%sv1=sqrt(vv1);sv2=sqrt(vv2);

%rectangle('position',[m1-sv1 m2-sv2 2*sv1 2*sv2]);
%axis([m1-(5*sv1) m1+(5*sv1) m2-(5*sv2) m2+(5*sv2)]);
%hold on
%for i=2:length(KG)
 %   if(io_const(i)==1)
  %      hold on
   %     plot(c_ratio(i),io_ratio(i),'ro', 'MarkerSize',2);
    %    hold on
    %end
    %if(io_const(i)==0)
     %   hold on
      %  plot(c_ratio(i),io_ratio(i),'ko', 'MarkerSize',2);
       % hold on
    %end
%end


%%





