load('data_2023_04_06_18_01_5C0_velocitiesfields_no_mean.mat')
%%
savedata=1;
if savedata==1
    Filename = sprintf('data_%s_5C0_enstrophy_kinetic_energy.mat', datetime('now','Format','yyyy_MM_dd_HH_mm'));
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
% qrangelog=linspace(log10(2*pi/480),log10(2*pi/8),qnum);
% qrange=10.^qrangelog;
cons=s1*s2;
fr_rang=15:597;%[50:50:250,295,350:50:550,595];
% fr_rang=[20,50,200,300,597];%[50:50:250,295,350:50:550,595];
%%
%r is the positive of the bacteria and R is the radius of the rings
for fr_n=fr_rang
    disp(fr_n);
    Stats=VelocitiesFields.(sprintf('frame_%d', fr_n));
    Vf=table2array(Stats(:,5))*0.4;
    XY=table2array(Stats(:,2));
    XY=XY-[cut_c_edge,cut_r_edge]-[1,1];
    XY=XY*0.4;
    N=size(XY,1);
    [qx,qy]=meshgrid(qrangex,qrangey);
    q2=qx.^2+qy.^2;
    q1=sqrt(q2);
    Eq_en=zeros(1,qnum);
    KEq=Eq_en;
    Fvx=zeros(qlengthy,qlengthx);
    Fvy=Fvx;
    tic
    for l = 1:N
        eqr=exp(-1i*(qx*XY(l,1)+qy*XY(l,2)));
        Fvx=Fvx+eqr*Vf(l,1);
        Fvy=Fvy+eqr*Vf(l,2);
    end
    Fvx=Fvx/N;
    Fvy=Fvy/N;
    toc
    cnt=0;
    for q=qrange
        cnt=cnt+1;
        Wq=8/3*pi*q*(3*hx);
        Wq2=4*(3*hx)/3;
        wq=1-(q1-q).^2/(3*hx)^2;
        wq(wq<0)=0;
        vq=Fvx.*qy-Fvy.*qx;
        Eq_en(cnt)=cons/Wq*sum((vq.*conj(vq)).*wq,'all');
        KEq(cnt)=sum((Fvx.*conj(Fvx)+Fvy.*conj(Fvy)).*wq,'all')/Wq2;
    end
    Eq_ens.(sprintf('frame_%d', fr_n))=Eq_en;
    KEqs.(sprintf('frame_%d', fr_n))=KEq;
end
if savedata==1
    save(Filename,'Eq_ens','KEqs','qrange')
end
%%
load('data_2023_09_30_21_17_5C0_enstrophy_kinetic_energy.mat')
fr_rang=15:597;
%% Figure for single frames

figkin=figure;
% Set units to inches
figkin.PaperUnits = 'inches';
% Set figure size
figkin.PaperPosition = [0 0 6 6];  % [left bottom width height]
hold on
col = jet(8);
cnt=0;
for fr_n=fr_rang([50,70,100,200,350,570]-14)%([25,45,55,70,100,200,350,570]-14)
    cnt=cnt+1;
    KEq=KEqs.(sprintf('frame_%d', fr_n));
    plot(log10(qrange(1:137)),log10(KEq(1:137)),'LineWidth',2,'Color', col(cnt,:))
    %writematrix([qrange(1:137); KEq(1:137)],'Kinetic_Energy.xlsx','Sheet',cnt)
    %     pt0=[x0,y0];
    %     pt1=[x0,y1];
    %     arr=pt1-pt0;
    %     quiver(pt0(1),pt0(2),arr(1),arr(2),'linewidth',2,'AutoScale','off','MaxHeadSize',1)
    %     x0=log10(qprange(arrowposition+cnt));
    %     y0=log10(Eq_el(arrowposition+cnt));
    %     plot(qprange,log10(Eq_el),'LineWidth',2)
end
legend([num2str(fr_rang(25-14)/5)],[num2str(fr_rang(70-14)/5)],[num2str(fr_rang(100-14)/5)], ...
   [num2str(fr_rang(200-14)/5)],[num2str(fr_rang(350-14)/5)],[num2str(fr_rang(570-14)/5)],'Location','southeast')
% legend([num2str(fr_rang(25-14)/5)], ...
%     [num2str(fr_rang(45-14)/5)],[num2str(fr_rang(55-14)/5)],[num2str(fr_rang(70-14)/5)],[num2str(fr_rang(100-14)/5)],[num2str(fr_rang(200-14)/5)],[num2str(fr_rang(350-14)/5)],[num2str(fr_rang(570-14)/5)],'Location','southeast')
 hold off
xlow=-1.6;
xhi=-0.5;
xlim([xlow,xhi])
xt=xlow:0.1:xhi;
xticks(xt)
xtlabel=cell(size(xt));
for i=1:length(xt)
    %         xlable{1,i}=['10^{',num2str(xt(i)),'}'];
    xtlabel{1,i}=['10^{',num2str(xt(i)),'}\newline',num2str(round(2*pi/10^(xt(i)),0))];
end
xticklabels(xtlabel)
ylow=-0.5;
yhi=2.5;
ylim([ylow,yhi])
yt=ylow:0.5:yhi;
yticks(yt)
ytlabel=cell(size(yt));
for i=1:length(yt)
    ytlabel{1,i}=['10^{',num2str(yt(i)),'}'];
end
yticklabels(ytlabel)
% set(findall(gcf,'-property','FontSize'),'FontSize',10,  'linewidth',2, 'fontweight', 'normal' )
% axis equal
% grid on
xlabel('q in (\mum^{-1}) \newline r in (\mum)','FontName','Times New Roman','FontSize',15)
ylabel('Kinetic Energy |K(q)|','FontName','Times New Roman','FontSize',24)
legend('Location', 'best');
legend('Orientation','vertical');
leg = legend('show');
title(leg,'Time (s)')
% set(gcf, 'Position',  [400, 400, 1000, 1000])
% ax = gca;  % Get the current axes handle
% ax.GridLineStyle = '--';  % Set grid line style to dotted
% ax.GridColor = [0 0 0];  % Set grid color to black
% set(0, 'DefaultAxesGridLineStyle', '--');
% set(0, 'DefaultAxesGridColor', [0 0 0]);
saveas(gcf,'figkin.png','png');
%==============================
%% title('Enstrophy Spectrum')
figenst=figure;
hold on
cnt=0;
for fr_n=fr_rang([25,70,100,200,350,570]-14)
    cnt=cnt+1;
    Eq_en=Eq_ens.(sprintf('frame_%d', fr_n));
    Eq_en(isnan(Eq_en))=0;
    plot(log10(qrange(1:137)),log10(Eq_en(1:137)),'LineWidth',2,'Color', col(cnt,:))
end
legend([num2str(fr_rang(25-14)/5)],[num2str(fr_rang(70-14)/5)],[num2str(fr_rang(100-14)/5)], ...
    [num2str(fr_rang(200-14)/5)],[num2str(fr_rang(350-14)/5)],[num2str(fr_rang(570-14)/5)],'Location','southeast')
hold off
xlow=-1.6;
xhi=-0.5;
xlim([xlow,xhi])
xt=xlow:0.1:xhi;
xticks(xt)
xtlabel=cell(size(xt));
for i=1:length(xt)
    xtlabel{1,i}=['10^{',num2str(xt(i)),'}\newline',num2str(round(2*pi/10^(xt(i)),0)),' \mum'];
end
xticklabels(xtlabel)
ylow=2;
yhi=6;
ylim([ylow,yhi])
yt=ylow:0.5:yhi;
yticks(yt)
ytlabel=cell(size(yt));
for i=1:length(yt)
    ytlabel{1,i}=['10^{',num2str(yt(i)),'}'];
end
yticklabels(ytlabel)
ylabel('Enstrophy Spectrum |{\Omega}(q)|','FontName','Times New Roman')
grid on
xlabel('q in (\mum^{-1}) \newline r in (\mum)')
set(findall(gcf,'-property','FontSize'),'FontSize',10)
legend('Location', 'northoutside');
legend('Orientation','vertical');
leg = legend('show');
title(leg,'Time (s)')
% set(gcf, 'Position',  [400, 400, 1000, 1000])
saveas(gcf,'figenst.png','png');
%% average
fr_ave=5;
col = jet(8);
n_ave=fr_ave+1;%number of frames used for averaging
fr_sel=fr_rang([25,45,55,70,100,200,350,570]-14);
for fr_n=fr_sel
    Eq_en_ave=0;
    for i = fr_n-fr_ave:fr_n+fr_ave
        temp=Eq_ens.(sprintf('frame_%d', i));
        temp(isnan(temp))=0;
        Eq_en_ave=Eq_en_ave+temp;
    end
    Eq_en_aves.(sprintf('frame_%d', fr_n))=Eq_en_ave/n_ave;
end
avenst=figure;
hold on
cnt=0;
for fr_n=fr_sel
    cnt=cnt+1;
    Eq_ave_en=Eq_en_aves.(sprintf('frame_%d', fr_n));
    plot(log10(qrange(1:137)),log10(Eq_ave_en(1:137)),'LineWidth',2,'Color', col(cnt,:))
    writematrix([qrange(1:137); Eq_ave_en(1:137)],'Average_Enstrophy_Spectrum.xlsx','Sheet',cnt)
end
% legend(num2str(fr_sel(1)/5),num2str(fr_sel(2)/5),num2str(fr_sel(3)/5),num2str(fr_sel(4)/5), ...
%    num2str(fr_sel(5)/5),num2str(fr_sel(6)/5))
% hold off
legend([num2str(fr_rang(25-14)/5)], ...
    [num2str(fr_rang(45-14)/5)],[num2str(fr_rang(55-14)/5)],[num2str(fr_rang(70-14)/5)],[num2str(fr_rang(100-14)/5)],[num2str(fr_rang(200-14)/5)],[num2str(fr_rang(350-14)/5)],[num2str(fr_rang(570-14)/5)],'Location','southeast')
 hold off
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

ylow=3;
yhi=7.5;
ylim([ylow,yhi])
yt=ylow:0.5:yhi;
yticks(yt)
ytlabel=cell(size(yt));
for i=1:length(yt)
    ytlabel{1,i}=['10^{',num2str(yt(i)),'}'];
end
yticklabels(ytlabel)
% axis equal
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',11)
xlabel('q in \mum^{-1} \newline r in \mum','FontName','Times New Roman','FontSize',18)
ylabel('Average Enstrophy Spectrum $|\overline{\Omega}(q)|$','FontName','Times New Roman','FontSize',18,'interpreter','latex')
legend('Location', 'best');
legend('Orientation','vertical');
leg = legend('show');
title(leg,'Time (s)')
% set(gcf, 'Position',  [400, 400, 1000, 1000])
saveas(gcf,'avenst.svg','svg');
%%
fr_ave=5;
col = jet(10);
fr_sel=fr_rang([25,35,45,55,70,100,200,350,570]-14);
for fr_n=fr_sel
    KEq_ave=0;
    for i = fr_n-fr_ave:fr_n+fr_ave
        KEq_ave=KEq_ave+KEqs.(sprintf('frame_%d', i));
    end
    KEq_aves.(sprintf('frame_%d', fr_n))=KEq_ave/n_ave;
end
avekin=figure;
% title('Average Kinetic Energy')
hold on
cnt=1;
for fr_n=fr_sel
    cnt=cnt+1;
    KEq_ave=KEq_aves.(sprintf('frame_%d', fr_n));
    plot(log10(qrange(1:137)),log10(KEq_ave(1:137)),'LineWidth',2,'Color', col(cnt,:))
    writematrix([qrange(1:137); KEq_ave(1:137)],'Average_Kinetic_Energy.xlsx','Sheet',cnt)
end
legend([num2str(fr_rang(25-14)/5)],[num2str(fr_rang(35-14)/5)], ...
    [num2str(fr_rang(45-14)/5)],[num2str(fr_rang(55-14)/5)],[num2str(fr_rang(70-14)/5)],[num2str(fr_rang(100-14)/5)],[num2str(fr_rang(200-14)/5)],[num2str(fr_rang(350-14)/5)],[num2str(fr_rang(570-14)/5)],'Location','southeast')
 hold off
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
set(findall(gcf,'-property','FontSize'),'FontSize',11)
ylow=0;
yhi=3;
ylim([ylow,yhi])
yt=ylow:0.5:yhi;
yticks(yt)
ytlabel=cell(size(yt));
for i=1:length(yt)
    ytlabel{1,i}=['10^{',num2str(yt(i)),'}'];
end
yticklabels(ytlabel)
ylabel('Average Kinetic Energy $|\overline{K}(q)|$','FontName','Times New Roman','FontSize',20,'interpreter','latex')
xlabel('q in \mum^{-1} \newline r in \mum','FontName','Times New Roman','FontSize',20)
grid on
legend('Location', 'best');
legend('Orientation','vertical');
leg = legend('show');
title(leg,'Time (s)')
% set(gcf, 'Position',  [400, 400, 1000, 1000])
saveas(gcf,'avekin.svg','svg');
%% fr_sel=fr_rang([50,150,250,450,570]-14);
for fr_n=fr_sel
    KEq_ave=0;
    for i = fr_n-fr_ave:fr_n+fr_ave
        KEq_ave=KEq_ave+KEqs.(sprintf('frame_%d', i));
    end
    KEq_aves.(sprintf('frame_%d', fr_n))=KEq_ave/n_ave;
end
avekin=figure;
% title('Average Kinetic Energy')
hold on
cnt=1;
for fr_n=fr_sel
    cnt=cnt+1;
    KEq_ave=KEq_aves.(sprintf('frame_%d', fr_n));
  loglog((qrange(1:137)),(KEq_ave(1:137)),'LineWidth',2,'Color', col(cnt,:))
end
legend([num2str(fr_sel(1)/5)],[num2str(fr_sel(2)/5)],[num2str(fr_sel(3)/5)],[num2str(fr_sel(4)/5)], ...
    [num2str(fr_sel(5)/5)],'Location','northwest')
hold off
% xlow=-1.6;
% xhi=-0.5;
% xlim([xlow,xhi])
% xt=xlow:0.1:xhi;
% xticks(xt)
% xtlabel=cell(size(xt));
% for i=1:length(xt)
%     xtlabel{1,i}=['10^{',num2str(xt(i)),'}\newline',num2str(round(2*pi/10^(xt(i)),0))];
% end
% xticklabels(xtlabel)
% xlabel('q in \mum^{-1} \newline r in \mum')
% ylow=0;
% yhi=2.5;
% ylim([ylow,yhi])
% yt=ylow:0.5:yhi;
% yticks(yt)
% ytlabel=cell(size(yt));
% for i=1:length(yt)
%     ytlabel{1,i}=['10^{',num2str(yt(i)),'}'];
% end
% yticklabels(ytlabel)
ylabel('Average Kinetic Energy $|\overline{K}(q)|$','FontName','Times New Roman','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
grid on
legend('Location', 'northoutside');
legend('Orientation','horizontal');
leg = legend('show');
title(leg,'Time (s)')
% set(gcf, 'Position',  [400, 400, 1000, 1000])
saveas(gcf,'avekin.png','png');