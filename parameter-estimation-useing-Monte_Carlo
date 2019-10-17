%Monte-Carlo simulation experiment
% for establishing parameter values
% under which +ve steady state will be established
clc
clear all
%intial condition..............
x0=1*ones(1,8);% eigth is the number of variable
%time range..........
tspan=0:0.05:300;
l= length(tspan);
% random experiment...........................
no_exp=10;
% a samll positive number for establishing steady state (fixing values upto two decimal place)
eps = 0.001;
% setting the rang of ammonium inputs
s1_low=2;
s1_high=3;
s1=s1_low:0.5:s1_high;
%..............................................
%intial steady-state assignment
steady_state=zeros(l,8);
%...............................

% random inputs in Monte-carlo simulation
for j=1:length(s1)
    %...................................
    % Start of the Monte-Carlo Simulation Experiment
    for i=1:no_exp
        kcat = rand(1);
        lon = (1/s1(j))*kcat+rand(1);non=rand(1);kon=rand(1);
        mon = s1(j)*non+rand(1);lcat=rand(1);
    if(non/mon < (lon/kcat))% by nagumo condition
%.....................................
% system outputs
[T,X]=ode45(@(t,x) GDH_system(t,x,s1(j),lon,kon,non,lcat,mon,kcat),tspan,x0);
%............................................................................
        if(X > zeros(size(X)))%positiveness conditions
            for k=1:10 %steady state condition
                if(abs(X(l+1-k,:)-X(l-k,:))< eps*ones(size(X(l-k,:))))
                s1_p(i)=s1(j);lon_p(i)=lon; kon_p(i)=kon; non_p(i)=non;
                lcat_p(i)=lcat; mon_p(i)=mon; kcat_p(i)=kcat; steady_state(i,:)=X(l-1,:);
                else
                 s1_p(i)=s1(j);lon_p(i)=0; kon_p(i)=0; non_p(i)=0;
                 lcat_p(i)=0; mon_p(i)=0; kcat_p(i)=0;    
                end
            end
        else
        s1_p(i)=s1(j);lon_p(i)=0; kon_p(i)=0; non_p(i)=0;
        lcat_p(i)=0; mon_p(i)=0; kcat_p(i)=0;  
        end
    else
    s1_p(i)=s1(j);lon_p(i)=0; kon_p(i)=0; non_p(i)=0;
    lcat_p(i)=0; mon_p(i)=0; kcat_p(i)=0;
    end
    end
    %....................................................
    % printing the experiment results in a .xls file
    ii=0; % dumy index for checking
        for i=1:no_exp
            if(lon_p(i)>0) 
            jj=i;
            par(j,:)=[s1_p(j) lon_p(jj) kon_p(jj) non_p(jj) lcat_p(jj) mon_p(jj) kcat_p(jj) steady_state(jj,:)];
            end
           if(lon_p(i)==0)
            ii=ii+1;
            if(ii==no_exp)
                disp('random experiment number is not sufficient for the steady state paremeter estimation:')
                disp(s1_p(j))
                par(j,:)=[s1_p(j) lon_p(ii) kon_p(ii) non_p(ii) lcat_p(ii) mon_p(ii) kcat_p(ii) steady_state(ii,:)];
            end
          end
        end
   
end

xlswrite('GDH_system_par.xls',par);
