function [P,WSSE,kappa,A,amp,period]=undulation(stripe,local_width,complex,draw,figname)
% Inputs:
% stripe: a binary image, in the form of a matrix of size m*n. Waves in the stripe are in horizontal direction. Default size 200*1216.
% local_width: width of the interval for calculating the detected angles. Default value is 32 pixels.
% complex: 1 for using double angle in complex plane, 0 otherwise.
% draw: 0 for not producing any figure, 1 for produceing a figure of the input stripe with its best fitting wave.
% Outputs:
% P: 4 by 1 vecotor, whose elements are values of parameters m, theta_max, lambda and horizontal shift.
% WSSE: Optional. Weighted sum of squared error for best fitting.
% kappa and A are as in formular (18) and (19) of the paper UNDULATIONS OF LAMELLAR LIQUID CRYSTALS...
% amp: the amplitude of the best fitting wave
% period: period of the best fitting wave
% fig: Optional, it produces a figure if draw =1.
% Best fitting iterations will stop when the number of iterations exceeds 3000.

[width,leng]=size(stripe);
inter_n=floor(leng/local_width);
Para_local=zeros(inter_n,7);%calculating local order parameters
for i=1:inter_n
    stats_local=mass_center(1,stripe(:,(i-1)*local_width+1:i*local_width),0,0,1,1);
    [~,~,~,~,~,~,~,~,~,qS,tS,tS_n,as,n_mc]=scalar_order_para_3(stats_local,0);
    [~,op,ang]=tensor_proc(qS);
    [~,op1,ang1]=tensor_proc(tS/n_mc);
    [~,op2,ang2]=tensor_proc(tS_n/as);
    Para_local(i,1)=(i-1)*local_width+1;
    Para_local(i,2:7)=[op,ang,op1,ang1,op2,ang2];
end
% intervals for initial guess
TH=linspace(0,pi/4,20);
Diffj1=zeros(length(TH)^3,5);
Diffj2=Diffj1;
Diffj3=Diffj1;
cnt=0;
l = Para_local(:,1)+floor(local_width/2);
for m=linspace(0.001,1,10)
    m2=m^2;
    for thm=TH
        thm2=thm^2;
        kappa=thm2/2*(1+1/m2);
        lam_st=sqrt(2*kappa-thm^2)/0.2;
        for lam=linspace(lam_st,5*lam_st,10)
            A=sqrt((2*kappa-thm^2))/(2*lam);
            for hshift=0:10:100
                cnt=cnt+1;
                theta=thm*ellipj(A*(l-hshift),m)*180/pi;
                diff2=Para_local(:,2).*(theta-Para_local(:,3)).^2;
                diff2s=sum(diff2);%op as weights
                Diffj1(cnt,:)=[-diff2s,m,thm,lam,hshift];
                diff2=Para_local(:,4).*(theta-Para_local(:,5)).^2;
                diff2s=sum(diff2);%op as weights
                Diffj2(cnt,:)=[-diff2s,m,thm,lam,hshift];
                diff2=Para_local(:,6).*(theta-Para_local(:,7)).^2;
                diff2s=sum(diff2);%op as weights
                Diffj3(cnt,:)=[-diff2s,m,thm,lam,hshift];
            end
        end
    end
end
min_dj1=max(Diffj1(:,1));
min_diff_indj1=find(Diffj1(:,1)==min_dj1);
min_dj2=max(Diffj2(:,1));
min_diff_indj2=find(Diffj2(:,1)==min_dj2);
min_dj3=max(Diffj3(:,1));
min_diff_indj3=find(Diffj3(:,1)==min_dj3);
%
Min_dj=zeros(3,5);
Min_dj(1,:)=Diffj1(min_diff_indj1(1),:);
Min_dj(2,:)=Diffj2(min_diff_indj2(1),:);
Min_dj(3,:)=Diffj3(min_diff_indj3(1),:);
[~,minid]=max(Min_dj(:,1));
% Searching for local minimal
w=Para_local(:,2*minid);
th=Para_local(:,2*minid+1);
if complex==1
    fun2 = @(X)jacobellip2_comp(X,l,th,w);
else
    fun2 = @(X)jacobellip2(X,l,th,w);
end
m0=Min_dj(minid,2);
thm0=Min_dj(minid,3);
lam0=Min_dj(minid,4);
hshift0=Min_dj(minid,5);
[P,WSSE] = fmincon(fun2,[m0,thm0,lam0,hshift0]',[],[],[],[],[0 -5*pi/12 0 0]',[1 5*pi/12 2000 2000]'); %constrained minimization
L = 1:leng;
m=P(1);
thm=P(2);
lam=P(3);
hshift=P(4);
kappa=thm^2/2*(1+1/m^2);
A=sqrt(2*kappa-thm^2)/(2*lam);
[~,cn,dn]=ellipj(A*(L-hshift),m);
U = thm/(A*sqrt(m))*log(dn-sqrt(m)*cn);
% period
if U(2)-U(1)>0
    for i = 2:leng-1
        if U(i+1)-U(i)<0
            break
        end
    end
    mark1=i;
    for j=mark1:leng-1
        if U(j+1)-U(j)>0
            break
        end
    end
    mark2=j;
else
    for i = 2:leng-1
        if U(i+1)-U(i)>0
            break
        end
    end
    mark1=i;
    for j=mark1:leng-1
        if U(j+1)-U(j)<0
            break
        end
    end
    mark2=j;
end
period=2*(mark2-mark1);
%
Umax=max(U);
Umin=min(U);
amp=(Umax-Umin)/2;
% figure
if draw == 1
    if minid==1
        str=('Weighted');
    else
        if minid==2
            str=('Unweighted');
        else
            str=('Weighted by area');
        end
    end
    const=amp-Umax;
    u=zeros(width,leng);
    for i = 1:leng
        j = round(U(i),0)+round(const,0)+ceil(width/2);
        u(j,i)=1;
    end
    u=flip(u);
    
    fig =figure('visible','off');
    cc=imfuse(stripe,u);
    imshow(cc);
    axis ij
    xticks(1:100:leng)
    xticklabels(num2cell((1:100:leng)-1))
    yticks([1:50:width,width])    
    yticklabels(num2cell([1:50:width,width]-width/2-1))
    axis on
    str_ver=['Stripe Width = ', num2str(width),newline,...
        'Amplitute = ',num2str(round(amp,2)),', Period = ',num2str(round(period,2)),newline,...
        'm = ',num2str(round(m,2)),', \theta_{max} = ' num2str(round(thm,2)),', \lambda = ' num2str(round(lam,2)),newline,...
        'Phase shift = ',num2str(round(hshift,2)),', Error = ',num2str(round(WSSE,1)),newline,...
        'Order Parameter Type = ',str];
    text(0,-0.75,str_ver,'FontSize',14, 'Units', 'Normalized', 'VerticalAlignment', 'Bottom')
    set(fig,'Position',[1,317,1440,480])
    saveas(fig,figname);
end
end
% main function ends
function [stats_t_r,stats_t_inn,Stats_r,Stats_inn,n_r,n_c,n_r_r,n_c_r,n_r_inn,n_c_inn,n_center_r,n_center_inn]=mass_center(fr_n,video,cut_r_edge,cut_c_edge,st_row,st_col)
imgb=logical(video(:,:,fr_n));
% imgb=imgb(cut_row_edge+1:end-cut_row_edge,cut_col_edge+1:end-cut_col_edge);
%imgb(imgb>0)=1;
%
% imgb = bwareaopen(imgb,6);%Remove all object containing fewer than 6 pixels.
[n_r,n_c] = size(imgb);
n_r_r = n_r-2*cut_r_edge;
n_c_r = n_c-2*cut_c_edge;
n_r_inn=n_r_r-2*(st_row-1);
n_c_inn=n_c_r-2*(st_col-1);
%CC = bwconncomp(imgb);
%
stats_t_r = regionprops('table',imgb,'Centroid','Orientation','PixelList','Area','MajorAxisLength','MinorAxisLength');
Stats_r = cat(2,stats_t_r.Centroid,stats_t_r.Orientation,stats_t_r.Area,stats_t_r.MajorAxisLength,stats_t_r.MinorAxisLength);
Stats_r(:,1:2) = Stats_r(:,1:2)-[cut_c_edge cut_r_edge];
[r,~]=find(Stats_r(:,1:2)<1);
Stats_r(r,:)=[];
stats_t_r(r,:)=[];
[r,~]=find(Stats_r(:,1)>n_c-2*cut_c_edge);
Stats_r(r,:)=[];
stats_t_r(r,:)=[];
[r,~]=find(Stats_r(:,2)>n_r-2*cut_r_edge);
Stats_r(r,:)=[];
stats_t_r(r,:)=[];
n_center_r=size(Stats_r,1);
%Find mc in the inner part.
Stats_inn = Stats_r;
stats_t_inn=stats_t_r;
[r,~]=find(Stats_inn(:,1)>st_col+n_c_inn-1);
Stats_inn(r,:)=[];
stats_t_inn(r,:)=[];
[r,~]=find(Stats_inn(:,1)<st_col);
Stats_inn(r,:)=[];
stats_t_inn(r,:)=[];
[r,~]=find(Stats_inn(:,2)>st_row+n_r_inn-1);
Stats_inn(r,:)=[];
stats_t_inn(r,:)=[];
[r,~]=find(Stats_inn(:,2)<st_row);
Stats_inn(r,:)=[];
stats_t_inn(r,:)=[];
n_center_inn=size(Stats_inn,1);
end

function [director,ave_angle,order_para,director1,ave_angle1,order_para1,director2,ave_angle2,order_para2,QS,TS,TS_n,A,i,Q,T,T_n]=scalar_order_para_3(stats,W) %scalar order parameter. Order_para1 is uniformly weighted.
Q=zeros(2,2); %Weighted
T=Q; %Unweighted
T_n=Q; %Weighted by number of pixels
A=0;%Total bact area
ratio=stats.MajorAxisLength./stats.MinorAxisLength;
for i = 1:size(stats,1) %the i-th mass center
    if ratio(i)>1
    x=stats.Centroid(i,1);
    y=stats.Centroid(i,2);
    Rx = stats.PixelList{i,1}(:,1)-x;% relative cooridinates of components of i
    Ry = y-stats.PixelList{i,1}(:,2);
    Qij = zeros(2,2);
    a=stats.Area(i);
    for j = 1:a
        Qij = Qij + [Rx(j),Ry(j)]'*[Rx(j),Ry(j)];
    end
    if W==0
        Q=Q+Qij;
    else
        xd=round(x,0);
        yd=round(y,0);
        w=W(xd,yd);
        Q=Q+w*Qij;
    end
    Qij_t=Qij/trace(Qij)-eye(2,2)/2;
    T=T+Qij_t;
    T_n=T_n+Qij_t*a;
    A=A+a;
    end
end

QS=Q;

[director,order_para,ave_angle]=tensor_proc(Q);

TS=T;
T=T/i; %num_list=trace(T)
[director1,order_para1,ave_angle1]=tensor_proc(T);

TS_n=T_n;
T_n=T_n/A; %num_list=trace(T)
[director2,order_para2,ave_angle2]=tensor_proc(T_n);
end

function [dir,op,angle,V,L]=tensor_proc(Q)
tr=trace(Q);
if tr >= 10^-12
    Q=Q/tr-eye(2,2)/2;
end
[V,L] = eig(Q);
dir=V(:,2);
op=2*L(2,2);
angle=atan(dir(2)/dir(1))*180/pi;
end

function f = jacobellip2(X,l,th,w)%X=[m,theta_max,lambda,hshift]
m=X(1);
m2=m^2;
thm=X(2);
thm2=thm^2;
lam=X(3);
k=thm2/2*(1+1/m2);
hshift=X(4);
A=sqrt(2*k-thm^2)/(2*lam);
theta=thm*ellipj(A*(l-hshift),m)*180/pi;
D2=w.*(th-theta).^2;
f=sum(D2,'all');
end

function f = jacobellip2_comp(X,l,th,w)%X=[m,theta_max,lambda,hshift]
th=th*pi/180;
m=X(1);
m2=m^2;
thm=X(2);
thm2=thm^2;
lam=X(3);
k=thm2/2*(1+1/m2);
hshift=X(4);
A=sqrt(2*k-thm^2)/(2*lam);
theta=thm*ellipj(A*(l-hshift),m);
D2=exp(1i*2*th)-exp(1i*2*theta);
D2=vecnorm(D2,2,2);
D2=w.*D2;
f=sum(D2,'all');
end