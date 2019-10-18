clear all;
clc;
mydata = xlsread('robust_matrix_random_exp1.xls');
x= mydata(:,4)./mydata(:,5);
c_r = mydata(:,11); p_r = mydata(:,12);
j=1;
for i=1:length(x)
    if(x(i)>1.2)
      x_s(j)=x(i);
      c_r_s(j)=c_r(i);
      p_r_s(j)=p_r(i);
      j=j+1;
    end
end
sub_rand_exp1_nadp_not_zero = [x_s' c_r_s' p_r_s'];
%..plot(f,cdate,pop)
mydata1 = xlsread('robust_matrix_robustfile1.xlsx');
.......................................................
   
subplot(1,2,1);
f=fit(x_s',c_r_s','poly1');
f1=fit(x_s',p_r_s','poly1');
plot(x,mydata(:,11),'bs');
niceplot;
hold on;
plot(f,'k--',x_s',c_r_s','bs');
niceplot;
hold on;
plot(x,mydata(:,12),'rs');
plot(f1,'k--',x_s',p_r_s','rs');
legend('off')
%b = gca; legend(b,'off');
%grid(ax1,'on')
%ax1.GridAlpha = 0.075;
niceplot;
subplot(1,2,2);
plot(mydata1(:,1),mydata1(:,6),'rs-');
niceplot;
hold on;
plot(mydata1(:,1), mydata1(:,5),'bs-');
niceplot;
hold off;
%legend('[GDH.NADP+.2OG]/[GDH.NADPH.2OG]','[Glu]/[NH4+]');
%legend({'[GDH.NADP^{+}.2OG]/[GDH.NADPH.2OG]','[Glu]/[NH4^{+}]'},'Location','northoutside');
%Box('off');
saveas(gcf, 'test123', 'epsc');
