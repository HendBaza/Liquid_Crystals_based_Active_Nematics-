% load('data_2023_01_30_tracking_5C0_all.mat')
% load('data_2023_04_20_16_05_5C0_vcorr_mean_subtracted.mat', 'VCorrws')
load('data_2023_04_06_18_01_5C0_velocitiesfields_no_mean.mat')
%%
importVF=VelocitiesFields_no_mean;
% t=0.2;%second between two frames
% T=2*t;
%1 pixel/2frames=0.4micron/0.4sec=1;
savedata=0;
ringwidth=16;
t2=ringwidth/2;
s1=1200;
s2=1800;
cut_r_edge=8;
cut_c_edge=12;
rmax=round(min(s1,s2)/4,0);
fr_rang=[93,95,97];%for 10C_0. % [50:50:250,295,350:50:550,595];%for 5C_0
%% <v(r),v(r+R)>
%r is the location of a bacterium and R is the radius of the rings. "< >"
%will average over (r), thus it becomes a function of R.
for fr_n=fr_rang%(2:end)
    disp(fr_n);
    Stats=importVF.(sprintf('frame_%d', fr_n));
    XY=table2array(Stats(:,1));
    XY=XY-[cut_c_edge,cut_r_edge]-[1,1];% adjust the center region such that (13,9), the left corner, is shifted to (0,0)
    V=table2array(Stats(:,2))*0.4;%in um/sec. 1 pixel/sec=0.4micron/sec;
    Vf=table2array(Stats(:,5))*0.4;
    N=size(XY,1);
    Distances=zeros(N);%in pixels
    Distances_shear=Distances;
    Distances_per=Distances;
    for i=1:N-1 %i-th bacterium in the frame
        xi=XY(i,1);
        yi=XY(i,2);
        for j=i+1:N
            xj=XY(j,1);
            yj=XY(j,2);%y-coordinates are up side down.
            v=[xj-xi,yj-yi];
            d=norm(v);
            Distances(i,j)=d;
            Distances(j,i)=d;
     
            if abs(v(2))<=ringwidth %only consider those whose distance in the shear direction to i is less than the width of the ring
                d=norm(v);
                Distances_per(i,j)=d;
                Distances_per(j,i)=d;
            end


            if abs(v(1))<=ringwidth
                d=norm(v);
                Distances_shear(i,j)=d;
                Distances_shear(j,i)=d;
            end
            
        end
    end

%     Nb=zeros(N,rmax);%number of bacteria in the ring centered at ith bacterium
%     Nbw=Nb;
%     Vdot=zeros(N,rmax);
%     Vdotn=Vdot;
%     Vdotw=Vdot;
%     Vdotnw=Vdot;
% 
%     Vfdot=Vdot;
%     Vfdotn=Vdot;
%     Vfdotw=Vdot;
%     Vfdotnw=Vdot;
%     
%     V2b=Nb;
%     Vf2b=Nb;
    %%
    tic
%     for i = 1:N
%         %         disp(i)
%         for r_in=1:rmax
%             r_mid=r_in+ringwidth/2;
%             r_out=r_in+ringwidth;
%             binring=Distances(:,i)>=r_in & Distances(:,i)<=r_out;
%             d=Distances(i,:);
%             d=d(binring);
%             d=abs(r_mid-d);
%             if isempty(d)~=1
%                 wk=1-4*d.^2/ringwidth^2;
%             else
%                 wk=0;
%             end
%             vdot=V(i,:)*V(binring,:)';
%             v2=V(i,:)*V(i,:)';
%             Vdot(i,r_in)=sum(vdot);%dot product of the vel of bacteria for i-th bacterium at r_in with i-th bacterium
%             Vdotn(i,r_in)=sum(vdot/v2);%normalized
%             Vdotw(i,r_in)=sum(vdot.*wk);
%             Vdotnw(i,r_in)=sum(vdot/v2.*wk);
% 
%             vfdot=Vf(i,:)*Vf(binring,:)';
%             vf2=Vf(i,:)*Vf(i,:)';
%             Vfdot(i,r_in)=sum(vfdot);%dot product of the vel of bacteria for i-th bacterium at r_in with i-th bacterium
%             Vfdotn(i,r_in)=sum(vfdot/vf2);%normalized
%             Vfdotw(i,r_in)=sum(vfdot.*wk);
%             Vfdotnw(i,r_in)=sum(vfdot/vf2.*wk);
% 
%             V2b(i,r_in)=norm(V(binring,:),'fro')^2;
%             Vf2b(i,r_in)=norm(Vf(binring,:),'fro')^2;
%             Nb(i,r_in)=sum(binring);%# of bacteria in r_in centered at i-th bacterium
%             Nbw(i,r_in)=sum(wk);
%         end
%     end
    [Nb,Nbw,Vdot,Vdotn,Vdotw,Vdotnw,Vfdot,Vfdotn,Vfdotw,Vfdotnw,V2b,Vf2b]=VCF(N,rmax,V,Vf,Distances,ringwidth);
    toc
    tic
    [Nb_shear,Nbw_shear,Vdot_shear,Vdotn_shear,Vdotw_shear,Vdotnw_shear,Vfdot_shear,Vfdotn_shear,Vfdotw_shear,Vfdotnw_shear,V2b_shear,Vf2b_shear]=VCF(N,rmax,V,Vf,Distances_shear,ringwidth);
    toc
    tic
    [Nb_per,Nbw_per,Vdot_per,Vdotn_per,Vdotw_per,Vdotnw_per,Vfdot_per,Vfdotn_per,Vfdotw_per,Vfdotnw_per,V2b_per,Vf2b_per]=VCF(N,rmax,V,Vf,Distances_per,ringwidth);
    toc
    Vdots.(sprintf('frame_%d', fr_n))=Vdot;
    Vdotws.(sprintf('frame_%d', fr_n))=Vdotw;
    Vdotns.(sprintf('frame_%d', fr_n))=Vdotn;
    Vdotnws.(sprintf('frame_%d', fr_n))=Vdotnw;
    Vfdots.(sprintf('frame_%d', fr_n))=Vfdot;
    Vfdotws.(sprintf('frame_%d', fr_n))=Vfdotw;
    Vfdotns.(sprintf('frame_%d', fr_n))=Vfdotn;
    Vfdotnws.(sprintf('frame_%d', fr_n))=Vfdotnw;
    Nbs.(sprintf('frame_%d', fr_n))=Nb;
    Nbws.(sprintf('frame_%d', fr_n))=Nbw;

    Vdots_shear.(sprintf('frame_%d', fr_n))=Vdot_shear;
    Vdotws_shear.(sprintf('frame_%d', fr_n))=Vdotw_shear;
    Vdotns_shear.(sprintf('frame_%d', fr_n))=Vdotn_shear;
    Vdotnws_shear.(sprintf('frame_%d', fr_n))=Vdotnw_shear;
    Vfdots_shear.(sprintf('frame_%d', fr_n))=Vfdot_shear;
    Vfdotws_shear.(sprintf('frame_%d', fr_n))=Vfdotw_shear;
    Vfdotns_shear.(sprintf('frame_%d', fr_n))=Vfdotn_shear;
    Vfdotnws_shear.(sprintf('frame_%d', fr_n))=Vfdotnw_shear;
    Nbs_shear.(sprintf('frame_%d', fr_n))=Nb_shear;
    Nbws_shear.(sprintf('frame_%d', fr_n))=Nbw_shear;

    Vdots_per.(sprintf('frame_%d', fr_n))=Vdot_per;
    Vdotws_per.(sprintf('frame_%d', fr_n))=Vdotw_per;
    Vdotns_per.(sprintf('frame_%d', fr_n))=Vdotn_per;
    Vdotnws_per.(sprintf('frame_%d', fr_n))=Vdotnw_per;
    Vfdots_per.(sprintf('frame_%d', fr_n))=Vfdot_per;
    Vfdotws_per.(sprintf('frame_%d', fr_n))=Vfdotw_per;
    Vfdotns_per.(sprintf('frame_%d', fr_n))=Vfdotn_per;
    Vfdotnws_per.(sprintf('frame_%d', fr_n))=Vfdotnw_per;
    Nbs_per.(sprintf('frame_%d', fr_n))=Nb_per;
    Nbws_per.(sprintf('frame_%d', fr_n))=Nbw_per;
    %%
    % end
    % for fr_n=fr_rang
    Ntotal=sum(Nb);%# of bacteria in r_in
    Ntotals.(sprintf('frame_%d', fr_n))=Ntotal;
    Ntotalw=sum(Nbw);%# of bacteria in r_in
    Ntotalws.(sprintf('frame_%d', fr_n))=Ntotalw;

    V2total=sum(V2b);
    V2totals.(sprintf('frame_%d', fr_n))=V2total;
    V2corr=V2total./Ntotal;
    V2Corrs.(sprintf('frame_%d', fr_n))=V2corr;
    VCorr=sum(Vdot)./Ntotal;
    VnCorr=sum(Vdotn)./Ntotal;
    VCorrs.(sprintf('frame_%d', fr_n))=VCorr;
    VnCorrs.(sprintf('frame_%d', fr_n))=VnCorr;
    VCorrw=sum(Vdotw)./Ntotalw;
    VnCorrw=sum(Vdotnw)./Ntotalw;
    VCorrws.(sprintf('frame_%d', fr_n))=VCorrw;
    VnCorrws.(sprintf('frame_%d', fr_n))=VnCorrw;
    
    Vf2total=sum(Vf2b);
    Vf2totals.(sprintf('frame_%d', fr_n))=Vf2total;
    Vf2corr=Vf2total./Ntotal;
    Vf2Corrs.(sprintf('frame_%d', fr_n))=Vf2corr;
    VfCorr=sum(Vfdot)./Ntotal;
    VfnCorr=sum(Vfdotn)./Ntotal;
    VfCorrs.(sprintf('frame_%d', fr_n))=VfCorr;
    VfnCorrs.(sprintf('frame_%d', fr_n))=VfnCorr;
    VfCorrw=sum(Vfdotw)./Ntotalw;
    VfnCorrw=sum(Vfdotnw)./Ntotalw;
    VfCorrws.(sprintf('frame_%d', fr_n))=VfCorrw;
    VfnCorrws.(sprintf('frame_%d', fr_n))=VfnCorrw;
    %
    Ntotal_shear=sum(Nb_shear);%# of bacteria in r_in
    Ntotals_shear.(sprintf('frame_%d', fr_n))=Ntotal_shear;
    Ntotalw_shear=sum(Nbw_shear);%# of bacteria in r_in
    Ntotalws_shear.(sprintf('frame_%d', fr_n))=Ntotalw_shear;

    V2total_shear=sum(V2b_shear);
    V2totals_shear.(sprintf('frame_%d', fr_n))=V2total_shear;
    V2corr_shear=V2total_shear./Ntotal_shear;
    V2Corrs_shear.(sprintf('frame_%d', fr_n))=V2corr_shear;
    VCorr_shear=sum(Vdot_shear)./Ntotal_shear;
    VnCorr_shear=sum(Vdotn_shear)./Ntotal_shear;
    VCorrs_shear.(sprintf('frame_%d', fr_n))=VCorr_shear;
    VnCorrs_shear.(sprintf('frame_%d', fr_n))=VnCorr_shear;
    VCorrw_shear=sum(Vdotw_shear)./Ntotalw_shear;
    VnCorrw_shear=sum(Vdotnw_shear)./Ntotalw_shear;
    VCorrws_shear.(sprintf('frame_%d', fr_n))=VCorrw_shear;
    VnCorrws_shear.(sprintf('frame_%d', fr_n))=VnCorrw_shear;
    
    Vf2total_shear=sum(Vf2b_shear);
    Vf2totals_shear.(sprintf('frame_%d', fr_n))=Vf2total_shear;
    Vf2corr_shear=Vf2total_shear./Ntotal_shear;
    Vf2Corrs_shear.(sprintf('frame_%d', fr_n))=Vf2corr_shear;
    VfCorr_shear=sum(Vfdot_shear)./Ntotal_shear;
    VfnCorr_shear=sum(Vfdotn_shear)./Ntotal_shear;
    VfCorrs_shear.(sprintf('frame_%d', fr_n))=VfCorr_shear;
    VfnCorrs_shear.(sprintf('frame_%d', fr_n))=VfnCorr_shear;
    VfCorrw_shear=sum(Vfdotw_shear)./Ntotalw_shear;
    VfnCorrw_shear=sum(Vfdotnw_shear)./Ntotalw_shear;
    VfCorrws_shear.(sprintf('frame_%d', fr_n))=VfCorrw_shear;
    VfnCorrws_shear.(sprintf('frame_%d', fr_n))=VfnCorrw_shear;
    %
    Ntotal_per=sum(Nb_per);%# of bacteria in r_in
    Ntotals_per.(sprintf('frame_%d', fr_n))=Ntotal_per;
    Ntotalw_per=sum(Nbw_per);%# of bacteria in r_in
    Ntotalws_per.(sprintf('frame_%d', fr_n))=Ntotalw_per;

    V2total_per=sum(V2b_per);
    V2totals_per.(sprintf('frame_%d', fr_n))=V2total_per;
    V2corr_per=V2total_per./Ntotal_per;
    V2Corrs_per.(sprintf('frame_%d', fr_n))=V2corr_per;
    VCorr_per=sum(Vdot_per)./Ntotal_per;
    VnCorr_per=sum(Vdotn_per)./Ntotal_per;
    VCorrs_per.(sprintf('frame_%d', fr_n))=VCorr_per;
    VnCorrs_per.(sprintf('frame_%d', fr_n))=VnCorr_per;
    VCorrw_per=sum(Vdotw_per)./Ntotalw_per;
    VnCorrw_per=sum(Vdotnw_per)./Ntotalw_per;
    VCorrws_per.(sprintf('frame_%d', fr_n))=VCorrw_per;
    VnCorrws_per.(sprintf('frame_%d', fr_n))=VnCorrw_per;
    
    Vf2total_per=sum(Vf2b_per);
    Vf2totals_per.(sprintf('frame_%d', fr_n))=Vf2total_per;
    Vf2corr_per=Vf2total_per./Ntotal_per;
    Vf2Corrs_per.(sprintf('frame_%d', fr_n))=Vf2corr_per;
    VfCorr_per=sum(Vfdot_per)./Ntotal_per;
    VfnCorr_per=sum(Vfdotn_per)./Ntotal_per;
    VfCorrs_per.(sprintf('frame_%d', fr_n))=VfCorr_per;
    VfnCorrs_per.(sprintf('frame_%d', fr_n))=VfnCorr_per;
    VfCorrw_per=sum(Vfdotw_per)./Ntotalw_per;
    VfnCorrw_per=sum(Vfdotnw_per)./Ntotalw_per;
    VfCorrws_per.(sprintf('frame_%d', fr_n))=VfCorrw_per;
    VfnCorrws_per.(sprintf('frame_%d', fr_n))=VfnCorrw_per;
end
if savedata==1
    Filename = sprintf('data_%s_10C0_vcorr_mean_subtracted_2pi.mat', datestr(now,'yyyy_mm_dd_HH_MM'));
    save(Filename,'Ntotals','Ntotalws','Nbs','Nbws','V2totals','V2Corrs','VCorrs','VnCorrs','VCorrws','VnCorrws','Vdots','Vdotws','Vdotns','Vdotnws', ...
        'Vf2totals','Vf2Corrs','VfCorrs','VfnCorrs','VfCorrws','VfnCorrws','Vfdots','Vfdotws','Vfdotns','Vfdotnws')
    Filename = sprintf('data_%s_10C0_vcorr_mean_subtracted_shear.mat', datestr(now,'yyyy_mm_dd_HH_MM'));
    save(Filename,'Ntotals_shear','Ntotalws_shear','Nbs_shear','Nbws_shear','V2totals_shear','V2Corrs_shear','VCorrs_shear','VnCorrs_shear','VCorrws_shear','VnCorrws_shear','Vdots_shear','Vdotws_shear','Vdotns_shear','Vdotnws_shear', ...
        'Vf2totals_shear','Vf2Corrs_shear','VfCorrs_shear','VfnCorrs_shear','VfCorrws_shear','VfnCorrws_shear','Vfdots_shear','Vfdotws_shear','Vfdotns_shear','Vfdotnws_shear')
    Filename = sprintf('data_%s_10C0_vcorr_mean_subtracted_per.mat', datestr(now,'yyyy_mm_dd_HH_MM'));
    save(Filename,'Ntotals_per','Ntotalws_per','Nbs_per','Nbws_per','V2totals_per','V2Corrs_per','VCorrs_per','VnCorrs_per','VCorrws_per','VnCorrws_per','Vdots_per','Vdotws_per','Vdotns_per','Vdotnws_per', ...
        'Vf2totals_per','Vf2Corrs_per','VfCorrs_per','VfnCorrs_per','VfCorrws_per','VfnCorrws_per','Vfdots_per','Vfdotws_per','Vfdotns_per','Vfdotnws_per')
end
%%
rmax=300;
figure
cnt=0;
%     x = linspace(0,(rmax+ringwidth)*0.4);
%     ym= x*0;
%     yd = -0.01+x*0;
%     yu = .01+x*0;
for fr_n=fr_rang%[50,100,150,200,250,295]%fr_rang
    cnt=cnt+1;
    VC=VfCorrws_per.(sprintf('frame_%d', fr_n));
    % VC=VfCorrws_shear.(sprintf('frame_%d', fr_n));
    % VC=VfCorrws.(sprintf('frame_%d', fr_n));
    subplot(3,4,cnt)
    plot(((1:300)+8)*0.4,VC,'LineWidth',2)
    hold on
    plot(((1:300)+8)*0.4,zeros(1,300),'LineWidth',1)
    xlabel('radius R in \mum')
%     ylabel('|v_f|^2 in \mum^2/s^2')
        ylabel('V_fCF')
        ylim([-0.3,1.3])
        xlim([0,(rmax+8)*0.4])
        grid on
    titlename=['Frame ',num2str(fr_n)];
    title(titlename)
end
%% comparison of frames
figure
hold on
for fr_n=fr_rang% [50,100,200,295,595] %for 5C_0
    VC=VfCorrws_per.(sprintf('frame_%d', fr_n));
    % VC=VfCorrws_shear.(sprintf('frame_%d', fr_n));
    % VC=VfCorrws.(sprintf('frame_%d', fr_n));
    plot(((1:300)+8)*0.4,VC,'LineWidth',2)
    xlabel('Distance ($\mu m$)','interpreter','latex')
%     ylabel('|v_f|^2 in \mum^2/s^2')
        ylabel('WVCF perpendicular to the shear direction ($\mu m^2/s^2$)','interpreter','latex')
        ylim([-0.3,1.3])
        xlim([0,(rmax+8)*0.4])
        grid on
end
plot(((1:300)+8)*0.4,zeros(1,300),'LineWidth',1,'Color','black')
hold off
legend('10 sec','20 sec','40 sec','60 sec','120 sec')
%% checking velocity in the rings
rs=16;
figure
imshow(importedvideo(:,:,fr_n),[0,3])
hold on
quiver(Stats.Centroid(:,1),Stats.Centroid(:,2),V(:,1),V(:,2),'LineWidth',2,'color','yellow','AutoScaleFactor',1);
axis ij
axis equal
axis on
for i=1:300:N
    binring=Distances(:,i)>=1 & Distances(:,i)<=16;
    C=XY(binring,:)+[cut_c_edge,cut_r_edge]+[1,1];
    viscircles(C,rs,'Color','white','LineWidth',1);
end
hold off
%% Energy Spectrum
savedata=0;
W=sum(1-4*(-8:8).^2/16^2)*0.4;
dr=0.4;%micron
dt=pi/180;%radian
% use qrange from 2023_04_25_energy_spectrum codes
% qrangelog=log10(1/s1):0.05:log10(1/(0.4*7.5));%half length of a bacterium = 0.4*7.5
% qrange=10.^qrangelog;%1/(sqrt(1199*1799)*0.4):0.01:1/(sqrt(1200*1800/6000)*0.4);
%corresponding to approximately the size of the region and to the bacterial velocity field spacing
E=zeros(size(qrange));
R=(9:308)*0.4;
figure
cnt_s=0;
for fr_n=fr_rang
    Ntotalw=Ntotalws.(sprintf('frame_%d', fr_n));
%     Ntotal=Ntotals.(sprintf('frame_%d', fr_n));
    cnt_s=cnt_s+1;
    VC=VfCorrws.(sprintf('frame_%d', fr_n));
    cnt=0;
    thrange=(-179:1:180)*pi/180;
    for q=qrange
        cnt=cnt+1;
        E(cnt)=q/(2*pi)*sum(exp(-1i*q*R).*VC.*Ntotalw.*R*W*dt);
%         E(cnt)=q/(2*pi)*sum(exp(-1i*q*R).*VC.*Ntotalw.*R*W*dt);
%         E(cnt)=q/(2*pi)*sum(exp(-1i*q*R).*VC.*Ntotal.*R*dt*dr);
    end
    Es.(sprintf('frame_%d', fr_n))=E;
    %
    subplot(4,3,cnt_s)
    plot(qrange,abs(E),'black','LineWidth',2)
    xlabel('Wave number q (\mum^{-1})')
    ylabel('Energy spectrum |E(q)|')
    grid on
    title('Frame ',num2str(fr_n))
end
%%
figure
hold on
cnt_s=0;
for fr_n=fr_rang
    cnt_s=cnt_s+1;
%     VCorrw=VCorrws.(sprintf('frame_%d', fr_n));
    cnt=0;
    E=Efbs.(sprintf('frame_%d', fr_n));%generated by energy_spectrum codes.
    %
%     subplot(4,3,cnt_s)
%     plot(log10(qrange),log10(abs(E).*qrange.^2),'LineWidth',2)
    plot(log10(qrange),log10(abs(E)),'black','LineWidth',2)
%     xlabel('Wave number q/q_{\it{l}}','FontName','Times New Roman')
        xlabel('Wave number q (\mum^{-1})')
    ylabel('Energy spectrum |E(q)|','FontName','Times New Roman')
    xt=-2.5:0.5:0;
    xticks(xt)
    xlable=cell(size(xt));
    for i=1:length(xt)
        xlable{1,i}=['10^{',num2str(xt(i)),'}\newline',num2str(round(2*pi/10^(xt(i)),0)),' \mum'];
    end
    xticklabels(xlable)
%     ylim([2,5])
%     yt=2:0.5:5;
%     yticks(yt)
%     ylable=cell(size(yt));
%     for i=1:length(yt)
%         ylable{1,i}=['10^{',num2str(yt(i)),'}'];
%     end
%     yticklabels(ylable)
    grid on
    title('Frame ',num2str(fr_n))
end
hold off
%%
fr_n=fr_rang(12);
E=Es.(sprintf('frame_%d', fr_n));
logq=log10(qrange*6/(2*pi));
loge=log10(abs(E));

plotsec=3;
XQ=zeros(12,2);
k1=1;%1;%1;%1;%1;%1;%1;%1;%2;
k2=2;%2;%2;%2;%3;% 2;%3;%2;%5;
x1=logq(k1);
y1=loge(k1);
x2=logq(k2);
y2=loge(k2);
XQ(1,1:2)=[x1,y1];
k=(y2-y1)/(x2-x1); qselect=(-90:90)./8; d=abs(qselect-k); q=qselect(d==min(d)); qn=q*8;
yn=y1+(x2-x1)*q;
XQ(2,1:2)=[x2,yn];
XQ(3,:)=[q,qn];
if plotsec>=2
    k1=6;%7;%8;%7;%4;%5;%8;%7;%7;%4;
    k2=7;%8;%10;%9;%5;%6;%10;%9;%8;%7;
    x1=logq(k1);
    y1=loge(k1);
    x2=logq(k2);
    y2=loge(k2);
    XQ(4,1:2)=[x1,y1];
    k=(y2-y1)/(x2-x1); d=abs(qselect-k); q=qselect(d==min(d)); qn=q*8; yn=y1+(x2-x1)*q;
    XQ(5,1:2)=[x2,yn];
    XQ(6,:)=[q,qn];
    if plotsec>=3
        k1=9;%8;%6;
        k2=10;%10;%7;
        x1=logq(k1);
        y1=loge(k1);
        x2=logq(k2);
        y2=loge(k2);
        XQ(7,1:2)=[x1,y1];
        k=(y2-y1)/(x2-x1); d=abs(qselect-k); q=qselect(d==min(d)); qn=q*8; yn=y1+(x2-x1)*q;
        XQ(8,1:2)=[x2,yn];
        XQ(9,:)=[q,qn];
        if plotsec==4
            k1=10;
            k2=11;
            x1=logq(k1);
            y1=loge(k1);
            x2=logq(k2);
            y2=loge(k2);
            XQ(10,1:2)=[x1,y1];
            k=(y2-y1)/(x2-x1); d=abs(qselect-k); q=qselect(d==min(d)); qn=q*8; yn=y1+(x2-x1)*q;
            XQ(11,1:2)=[x2,yn];
            XQ(12,:)=[q,qn];
        end
    end
end
qs.(sprintf('frame_%d', fr_n))=XQ;
%
figure
plot(logq,loge,'black','LineWidth',2)
hold on
plot(XQ(1:2,1),XQ(1:2,2)+0.15,'red','LineWidth',2)
lg1=['\alpha=',num2str(XQ(3,1))];
legend('E',lg1,'Location','northwest')
if plotsec>=2
    plot(XQ(4:5,1)-0.02,XQ(4:5,2)-0.1,'blue','LineWidth',2)
    legend('E',lg1,lg2,'Location','northwest')
    lg2=['\alpha=',num2str(XQ(6,1))];
    if plotsec>=3
        plot(XQ(7:8,1)-0.01,XQ(7:8,2)+0.15,'r--','LineWidth',2)
        lg3=['\alpha=',num2str(XQ(9,1))];
        legend('E',lg1,lg2,lg3,'Location','northwest')
        if plotsec==4
            plot(XQ(10:11,1)-0.03,XQ(10:11,2)+0.02,'b--','LineWidth',2)
            lg4=['\alpha=',num2str(XQ(12,1))];
            legend('E',lg1,lg2,lg3,lg4,'Location','northwest')
        end
    end
end
hold off
xlabel('Wave number q \mum^{-1}','FontName','Times New Roman')
% xlabel('Wave number q/q_{\it{l}}','FontName','Times New Roman')
ylabel('Energy spectrum |E(q)|','FontName','Times New Roman')
xlim([-3,-0.75])
ylim([-1,2.5])
xticks(-2.5:0.5:-0.5)
xticklabels({'10^{-2.5}','10^{-2}','10^{-1.5}','10^{-1}','10^{-0.5}'})
yticks(-1:0.5:3)
yticklabels({'10^{-1}','10^{-0.5}','10^{0}','10^{0.5}' ...
    ,'10^{1}','10^{1.5}','10^{2}','10^{2.5}'})
grid on
title('Frame ',num2str(fr_n))