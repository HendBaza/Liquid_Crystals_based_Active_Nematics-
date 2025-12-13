load('data_2023_01_30_tracking_5C0_all.mat')
            %%
savedata=1;
if savedata==1
    Filename = sprintf('data_%s_5C0_elastic_energy.mat', datetime('now','Format','yyyy_MM_dd_HH_mm'));
end
s1=(1200)*0.4;
s2=(1800)*0.4;
cut_r_edge=8;
cut_c_edge=12;
Kx=0:38;
qrangex=2*pi/s2*Kx;
Ky=0:26;
qrangey=2*pi/s1*Ky;
qrangex=[-qrangex(end:-1:2),qrangex(1:end)];
qrangey=[-qrangey(end:-1:2),qrangey(1:end)];
hx=qrangex(2)-qrangex(1);
hy=qrangey(2)-qrangey(1);
qlengthx=length(qrangex);
qlengthy=length(qrangey);
qnum=300;
qrange=linspace(3*2*pi/s2,2*pi/10,qnum);
cons=s1*s2/2;
fr_rang=15:597;%[50:50:250,295,350:50:550,595];
%% %r is the positive of the bacteria and R is the radius of the rings
for fr_n=fr_rang
    disp(fr_n);
    Stats=TRACKING.(sprintf('frame_%d', fr_n));
    thetas=table2array(Stats(:,3));
    S=sin(thetas);
    C=cos(thetas);
    S2=sin(2*thetas);
    Q11=2*C.^2-1;
    Q12=S2;
    Q22=2*S.^2-1;
    XY=table2array(Stats(:,2));
    XY=XY-[cut_c_edge,cut_r_edge]-[1,1];
    %     XY=round(XY,0)-[cut_c_edge,cut_r_edge]-[1,1];% adjust the center region such that (13,9), the left corner, is shifted to (0,0)
    XY=XY*0.4;
    N=size(XY,1);
    [qx,qy]=meshgrid(qrangex,qrangey);
    q2=qx.^2+qy.^2;
    q1=sqrt(q2);
    F11=zeros(qlengthy,qlengthx);
    F12=F11;
    F22=F11;
    Eq_el=zeros(1,qnum);
    tic
    for l = 1:N
        F11=F11+exp(-1i*(qx*XY(l,1)+qy*XY(l,2)))*Q11(l);
        F12=F12+exp(-1i*(qx*XY(l,1)+qy*XY(l,2)))*Q12(l);
        F22=F22+exp(-1i*(qx*XY(l,1)+qy*XY(l,2)))*Q22(l);
    end
    F11=F11/N;
    F12=F12/N;
    F22=F22/N;
    E_el=cons*sum(F11.*(conj(F11)+2*F12.*conj(F12)+F22.*conj(F22)).*q2,'all');
    toc
%     E_els.(sprintf('frame_%d', fr_n))=E_el;
    cnt=0;
    for q=qrange
        cnt=cnt+1;
        Wq=8/3*pi*q*(3*hx);
        wq=1-(q1-q).^2/(3*hx)^2;
        wq(wq<0)=0;
        Eq_el(cnt)=cons/Wq*sum((F11.*conj(F11)+2*F12.*conj(F12)+F22.*conj(F22)).*q2.*wq,'all');
    end
    Eq_els.(sprintf('frame_%d', fr_n))=Eq_el;
end
if savedata==1
    save(Filename,'Eq_els','qrange')
end
%%
load('data_2023_09_30_13_05_5C0_elastic_energy.mat')
fr_rang=15:597;
%% Figure for single frames
figure
col = jet(10);
hold on
cnt=0;
for fr_n=fr_rang([15,25,35,45,55,70,100,200,350,570]-14)%fr_rang([25,70,100,200,350,570]-14)%
    cnt=cnt+1;
    Eq_el=Eq_els.(sprintf('frame_%d', fr_n));
    Eq_el(isnan(Eq_el))=0;
%     qslope=qrange(50:137);
%     plot(log10(qslope),log10(150000*qslope.^1.8),log10(qslope),log10(400000*qslope.^1.3),'LineWidth',2,'Color', 'k','LineStyle','-.')
%     hold on;
    plot(log10(qrange(1:137)),log10(Eq_el(1:137)),'LineWidth',2,'Color', col(cnt,:))
%     writematrix(qrange(1:137), Eq_el(1:137),'Elasticenergy.xlsx','Sheet',cnt)
%     pt0=[x0,y0];
%     pt1=[x0,y1];
%     arr=pt1-pt0;
%     quiver(pt0(1),pt0(2),arr(1),arr(2),'linewidth',2,'AutoScale','off','MaxHeadSize',1)
%     x0=log10(qprange(arrowposition+cnt));
%     y0=log10(Eq_el(arrowposition+cnt));
%     plot(qprange,log10(Eq_el),'LineWidth',2)
end
legend([num2str(fr_rang(15-14)/5)],[num2str(fr_rang(25-14)/5)],[num2str(fr_rang(35-14)/5)], ...
    [num2str(fr_rang(45-14)/5)],[num2str(fr_rang(55-14)/5)],[num2str(fr_rang(70-14)/5)],[num2str(fr_rang(100-14)/5)],[num2str(fr_rang(200-14)/5)],[num2str(fr_rang(350-14)/5)],[num2str(fr_rang(570-14)/5)],'Location','southeast')
hold off
grid on;
xlow=-1.6;
xhi=-0.5;
xlim([xlow,xhi])
xt=xlow:0.1:xhi;
xticks(xt)
xtlabel=cell(size(xt));
for i=1:length(xt)
    xtlabel{1,i}=['10^{',num2str(xt(i)),'}\newline',num2str(round(2*pi/10^(xt(i)),0))];
end
xticklabels(xtlabel)
ylow=2.5;
yhi=5;
ylim([ylow,yhi])
yt=ylow:0.5:yhi;
yticks(yt)
ytlabel=cell(size(yt));
for i=1:length(yt)
    ytlabel{1,i}=['10^{',num2str(yt(i)),'}'];
end
yticklabels(ytlabel)
set(findall(gcf,'-property','FontSize'),'FontSize',11,  'linewidth',2)
ylabel('Elastic Energy |E_{els}(q)|','FontName','Times New Roman','FontSize',20)
xlabel('q in (\mum^{-1}) \newline r in (\mum)','FontName','Times New Roman','FontSize',20)
legend('Location', 'best');
legend('Orientation','vertical');
leg = legend('show');
title(leg,'Time (s)')
% ax = gca;  % Get the current axes handle
% ax.GridLineStyle = '--';  % Set grid line style to dotted
% ax.GridColor = [0 0 0];  % Set grid color to black
% set(0, 'DefaultAxesGridLineStyle', '-.');
% set(0, 'DefaultAxesGridColor', [0 0 0]);
saveas(gcf,'ElasticEnergy.svg','svg');
%% Figure for single frames
% figure
% col = jet(10);
% hold on
% cnt=0;
% for fr_n=fr_rang([15,25,35,45,55,70,100,200,350,570]-14)%fr_rang([25,70,100,200,350,570]-14)%
%     cnt=cnt+1;
%     Eq_el=Eq_els.(sprintf('frame_%d', fr_n));
%     Eq_el(isnan(Eq_el))=0;
% %     y1=log10(Eq_el(arrowposition-1+cnt));
%     loglog(log10(qrange(1:137)),log10(Eq_el(1:137)),'LineWidth',2,'Color', col(cnt,:))
% %     pt0=[x0,y0];
% %     pt1=[x0,y1];
% %     arr=pt1-pt0;
% %     quiver(pt0(1),pt0(2),arr(1),arr(2),'linewidth',2,'AutoScale','off','MaxHeadSize',1)
% %     x0=log10(qprange(arrowposition+cnt));
% %     y0=log10(Eq_el(arrowposition+cnt));
% %     plot(qprange,log10(Eq_el),'LineWidth',2)
% end
% legend([num2str(fr_rang(15-14)/5)],[num2str(fr_rang(25-14)/5)],[num2str(fr_rang(35-14)/5)], ...
%     [num2str(fr_rang(45-14)/5)],[num2str(fr_rang(55-14)/5)],[num2str(fr_rang(70-14)/5)],[num2str(fr_rang(100-14)/5)],[num2str(fr_rang(200-14)/5)],[num2str(fr_rang(350-14)/5)],[num2str(fr_rang(570-14)/5)],'Location','southeast')
% % legend(['Frame ',num2str(fr_rang(25-14))],['Frame ',num2str(fr_rang(70-14))],['Frame ',num2str(fr_rang(100-14))], ...
% %     ['Frame ',num2str(fr_rang(200-14))],['Frame ',num2str(fr_rang(350-14))],['Frame ',num2str(fr_rang(570-14))],'Location','southeast')
% hold off
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% % xlow=-1.6;
% %     xhi=-0.5;
% %     xlim([xlow,xhi])
% %     xt=xlow:0.1:xhi;
% %     xticks(xt)
% %     xtlabel=cell(size(xt));
% %     for i=1:length(xt)
% % %         xlable{1,i}=['10^{',num2str(xt(i)),'}'];
% %         xtlabel{1,i}=['10^{',num2str(xt(i)),'}\newline',num2str(round(2*pi/10^(xt(i)),0)),' \mum'];
% %     end
% %     xticklabels(xtlabel)
% % ylow=3;%;-7;
% % yhi=5;%-4;
% % ylim([ylow,yhi])
% % yt=ylow:0.5:yhi;
% % yticks(yt)
% % ytlabel=cell(size(yt));
% % for i=1:length(yt)
% %     ytlabel{1,i}=['10^{',num2str(yt(i)),'}'];
% % end
% % yticklabels(ytlabel)
% % ylabel('Elastic Energy |E_{els}(q)|','FontName','Times New Roman')
% % xlabel('q in (\mum^{-1}) \newline r in (\mum)')
% % % axis equal
% % legend('Location', 'northoutside');
% % legend('Orientation','horizontal');
% % leg = legend('show');
% % title(leg,'Time (s)')
