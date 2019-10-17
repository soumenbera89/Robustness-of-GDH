function dx = GDH_system(t,x,k1,NADPH,k4,k2,k3,KG,k6,k5,NH4,k8,k9d,NADP,k7,k10d,k10,k12,k9,k11,k8d,k7d)
% External input of GDH system................
%% NADPH, KG, NH4
%%Zero column matrix.....
dx=zeros(8,1);
%GDH system............
%%x(1)=[GDH.NADPH],x(2)=[GDH.NADPH.KG],x(3)=[GDH.NADPH.KG.NH4-GDH.GLU.NADP+],x(4)=[GDH.NADP+]
%%x(5)=[GDH.NADP+.KG],x(6)=[GDH.GLU],x(7)=[GDH],x(8)=[GLU]
dx(1)=k1*x(7)*NADPH+k4*x(2)-k2*x(1)-k3*KG*x(1);
dx(2)=k3*KG*x(1)+k6*x(3)-k4*x(2)-k5*NH4*x(2);
dx(3)=k5*NH4*x(2)+k8*x(4)*x(8)+k9d*NADP*x(6)-(k6+k7+k10d)*x(3);
dx(4)=k7*x(3)+k10*NADP*x(7)+k12*x(5)-k8*x(8)*x(4)-k9*x(4)-k11*KG*x(4);
dx(5)=k11*KG*x(4)-k12*x(5);
dx(6)=k8d*x(8)*x(7)+k10d*x(3)-k7d*x(6)-k9d*NADP*x(6);
dx(7)=k2*x(1)+k7d*x(6)+k9*x(4)-k8d*x(8)*x(7)-k10*NADP*x(7)-k1*NADPH*x(7);
dx(8)=k7d*x(6)-k8*x(8)*x(4)-k8d*x(8)*x(7);
return
