% load('Jan242020_10um_5C0_1_300.mat')
importedvideo=Jan262020_82;
st=90;%11:299%<-------------------------------------------------------------------------------------------------range of frames to be analyzed
nd=101;
st_import=82;% first frame of the imported video
seq_adj=st_import-1;
%%
savedata=1;
time=0.4;% 0.2 second between two frames. We calculate the velocity by the displacement of bacteria in three frames;
thre_area=15; %min area of one bacterium
thre_ar=1.5; %min aspect ratio of one bacterium
thre_cl=12; %ignore bacteria whose mass centers are too close to the edges.
thre_r=8;
%generate tracking table
for fr_N_ab=st+1:nd-1
    fprintf('Analyzing Frame: ')
    sprintf('%0.0f',fr_N_ab)
    %<-----------------------------------------------------------------------------------------------choose a video
    % if fr_N_ab<=300
    %     X=v1_300g(:,:,fr_N_ab);
    % else
    %     X=v1_301(:,:,fr_N_ab-130);
    % end
    X=importedvideo(:,:,fr_N_ab-seq_adj);
    tic
    stats = regionprops('table',logical(X),'Centroid','Orientation','Area','MajorAxisLength','MinorAxisLength','PixelList');
    stats(stats.Area<thre_area,:)=[];
    stats(stats.Centroid(:,1)<=thre_cl | stats.Centroid(:,1)>1824-thre_cl,:)=[];
    stats(stats.Centroid(:,2)<=thre_r | stats.Centroid(:,2)>1216-thre_r,:)=[];
    stats.AspectRatio=stats.MajorAxisLength./stats.MinorAxisLength;
    stats(stats.AspectRatio<thre_ar,:)=[];
    stats.Orientation=-stats.Orientation*pi/180;%ranging from -pi/2 to pi/2. postive angles are from positive x-axis to positive y-axis.
    n_bac=size(stats,1);
    stats.Centroid_p=zeros(n_bac,2);%previous frame
    stats.Centroid_x=zeros(n_bac,2);%next frame
    stats.displacement_p=zeros(n_bac,1);%previous frame
    stats.displacement_x=zeros(n_bac,1);%next frame
    toc
    %
    fr_p=fr_N_ab-1;
    %<-----------------------------------------------------------------------------------------------choose a video
    % if fr_p<=300
    %     X_p=v1_300g(:,:,fr_p);
    % else
    %     X_p=v1_301(:,:,fr_p-130);
    % end
    X_p=importedvideo(:,:,fr_p-seq_adj);
    Stats_p = regionprops('table',logical(X_p),'Centroid','Orientation','Area','MajorAxisLength','MinorAxisLength','PixelList');
    Stats_p(Stats_p.Centroid(:,1)<=thre_cl | Stats_p.Centroid(:,1)>1824-thre_cl,:)=[];
    Stats_p(Stats_p.Centroid(:,2)<=thre_r | Stats_p.Centroid(:,2)>1216-thre_r,:)=[];
    Stats_p(Stats_p.Area<thre_area,:)=[];
    Stats_p(Stats_p.MajorAxisLength./Stats_p.MinorAxisLength<thre_ar,:)=[];
    Stats_p.Orientation=-Stats_p.Orientation*pi/180;
    fr_x=fr_N_ab+1;
    %<-----------------------------------------------------------------------------------------------choose a video
    % if fr_x<=300
    %     X_x=v1_300g(:,:,fr_x);
    % else
    %     X_x=v1_301(:,:,fr_x-130);
    % end
    X_x=importedvideo(:,:,fr_x-seq_adj);
    Stats_x = regionprops('table',logical(X_x),'Centroid','Orientation','Area','MajorAxisLength','MinorAxisLength','PixelList');
    Stats_x(Stats_x.Centroid(:,1)<=thre_cl | Stats_x.Centroid(:,1)>1824-thre_cl,:)=[];
    Stats_x(Stats_x.Centroid(:,2)<=thre_r | Stats_x.Centroid(:,2)>1216-thre_r,:)=[];
    Stats_x(Stats_x.Area<thre_area,:)=[];
    Stats_x(Stats_x.MajorAxisLength./Stats_x.MinorAxisLength<thre_ar,:)=[];
    Stats_x.Orientation=-Stats_x.Orientation*pi/180;
    tic
    for i=1:n_bac
        ctr=stats.Centroid(i,:);
        %tracking the previous and the current
        v_p=ctr-Stats_p.Centroid;
        tracked=0;
        cnt_trk=0;
        d_p=v_p(:,1).^2+v_p(:,2).^2;
        while tracked==0 && cnt_trk<=5
            cnt_trk=cnt_trk+1;
            mindp=min(d_p);
            % check if there are overlap pixels
            pix1=stats.PixelList{i,1};
            p1=pix1(:,1)*10000+pix1(:,2);
            ctr_p=Stats_p.Centroid(d_p==mindp,:);
            pix2=Stats_p.PixelList(d_p==mindp,:);
            pix2=pix2{1,1};
            p2=pix2(:,1)*10000+pix2(:,2);
            if isempty(intersect(p1,p2))==0
                stats.Centroid_p(i,:)=ctr_p;
                stats.displacement_p(i)=sqrt(mindp);
                tracked=1;
            else
                d_p(d_p==mindp)=Inf;
            end
        end
        %tracking the current and the next
        v_x=Stats_x.Centroid-ctr;%make center (0,0)
        tracked=0;
        cnt_trk=0;
        d_x=v_x(:,1).^2+v_x(:,2).^2;
        while tracked==0 && cnt_trk<=5
            cnt_trk=cnt_trk+1;
            mindx=min(d_x);
            %         mindxid=find(d_x==mindx);
            % check if there are overlap pixels
            pix1=stats.PixelList{i,1};
            p1=pix1(:,1)*10000+pix1(:,2);
            ctr_x=Stats_x.Centroid(d_x==mindx,:);
            pix2=Stats_x.PixelList(d_x==mindx,:);
            pix2=pix2{1,1};
            p2=pix2(:,1)*10000+pix2(:,2);
            if isempty(intersect(p1,p2))==0
                stats.Centroid_x(i,:)=ctr_x;
                stats.displacement_x(i)=sqrt(mindx);
                tracked=1;
            else
                d_x(d_x==mindx)=Inf;
            end
        end
    end
    toc
    stats(:,6)=[];
%     STATS.(sprintf('frame_%d', fr_N_ab))=stats;
    %
%     stats(stats.displacement_x>thre_disp,:)=[];
%     stats(stats.displacement_p>thre_disp,:)=[];
    TRACKING.(sprintf('frame_%d', fr_N_ab))=stats;
end
%% velocity field
for fr_N_ab=st+1:nd-1
    tic
    Stats_crt=TRACKING.(sprintf('frame_%d', fr_N_ab));
    Stats_crt(Stats_crt.Centroid_p(:,1)==0 | Stats_crt.Centroid_x(:,1)==0,:)=[];
    loc=Stats_crt.Centroid;
    N=size(Stats_crt,1);
    Distances=zeros(N);
    for i=1:N-1
        xi=loc(i,1);
        yi=loc(i,2);
        for j=i+1:N
            xj=loc(j,1);
            yj=loc(j,2);
            vel=[xj-xi,yj-yi];
            d=norm(vel);
            Distances(i,j)=d;
            Distances(j,i)=d;
        end
    end
    toc
    %
    tic
    V=table;
    V.Centroid=Stats_crt.Centroid;
    V.DisplaceVel=(Stats_crt.Centroid_x-Stats_crt.Centroid_p)/time;
    V.Orientation=Stats_crt.Orientation;
    V.Flow=zeros(N,2);
    V.Swim=zeros(N,2);
    V.BactN=zeros(N,1);% Number of the bacteria in the radius of 19.2 microns (48 pixels).
    for i = 1:N
        stats_k=V(Distances(:,i)<=48,:);
        K=size(stats_k,1);
        M=zeros(2*K,2);
        V_per=zeros(2*K,1);
        for k = 1:K
            phi=stats_k.Orientation(k);
            n1=cos(phi);
            n2=sin(phi);
            n=[n1;n2];
            dx=stats_k.DisplaceVel(k,1);
            dy=stats_k.DisplaceVel(k,2);
            v=[dx;dy];
            V_per((k-1)*2+1:(k-1)*2+2)=v-(n'*v)*n;%in pixels/sec
            M((k-1)*2+1,1)=1-n1*n1;
            M((k-1)*2+1,2)=-n1*n2;
            M((k-1)*2+2,1)=-n2*n1;
            M((k-1)*2+2,2)=1-n2*n2;
        end
        v_f=M\V_per;
        V.Flow(i,:)=v_f;
        V.BactN(i)=K;
    end
    V.Swim=V.DisplaceVel-V.Flow;
    VelocityField.(sprintf('frame_%d', fr_N_ab))=V;
    toc
end
%%
if savedata==1
    Filename = sprintf('data_%s_10C0_tracking_vel.mat', datestr(now,'yyyy_mm_dd_HH_MM'));
    save(Filename,'VelocityField')
end
%% Draw Fields of velocities
mag1=1;
fr_N_ab=95;%<-----------------------------------------------------------------------------------------------choose a frame to show
V=VelocityField.(sprintf('frame_%d', fr_N_ab));
figure
imshow(importedvideo(:,:,fr_N_ab-seq_adj)*3+importedvideo(:,:,fr_N_ab-seq_adj+1),[0,5])
hold on
quiver(V.Centroid(:,1),V.Centroid(:,2),mag1*V.DisplaceVel(:,1),mag1*V.DisplaceVel(:,2),'LineWidth',2,'color','blue','AutoScale','off')%<--from previous frame to the next.
quiver(V.Centroid(:,1),V.Centroid(:,2),mag1*V.Flow(:,1),mag1*V.Flow(:,2),'LineWidth',2,'color','red','AutoScale','off')
quiver(V.Centroid(:,1),V.Centroid(:,2),mag1*V.Swim(:,1),mag1*V.Swim(:,2),'LineWidth',2,'color','yellow','AutoScale','off')
axis on
axis ij
hold off
str_title=['Velocity Field of Frame ',num2str(fr_N_ab)];
title(str_title)
legend('Total','Flow','Swim','TextColor','white','Color','black')
% Draw single velocity
mag1=1;
% mag2=10;
%%
figure
% imshow(v1_300g(:,:,fr_N_ab),[0,1])
imshow(importedvideo(:,:,fr_N_ab-seq_adj)*3+importedvideo(:,:,fr_N_ab-seq_adj+1),[0,5])
hold on
% quiver(V.Centroid(:,1),V.Centroid(:,2),time*V.DisplaceVel(:,1),time*V.DisplaceVel(:,2),'LineWidth',2,'color','yellow','AutoScale','off')
quiver(V.Centroid(:,1),V.Centroid(:,2),mag1*V.DisplaceVel(:,1),mag1*V.DisplaceVel(:,2),'LineWidth',2,'color','yellow','AutoScale','on')
% quiver(V.Centroid(:,1),V.Centroid(:,2),mag1*V.Flow(:,1),mag1*V.Flow(:,2),'LineWidth',2,'color','yellow','AutoScale','off')
% quiver(V.Centroid(:,1),V.Centroid(:,2),mag1*V.Swim(:,1),mag1*V.Swim(:,2),'LineWidth',2,'color','yellow','AutoScale','on')
% quiver(V.Centroid(:,1),V.Centroid(:,2),mag2*cos(V.Orientation(:)),mag2*sin(V.Orientation(:)),'LineWidth',2,'color','red','AutoScale','off','ShowArrowHead','off')
% quiver(V.Centroid(:,1),V.Centroid(:,2),-mag2*cos(V.Orientation(:)),-mag2*sin(V.Orientation(:)),'LineWidth',2,'color','red','AutoScale','off','ShowArrowHead','off')
axis on
axis ij
hold off
str_title=['Velocity Field of Frame ',num2str(fr_N_ab)];
title(str_title)
%% trajectory visual checking
%from current to the next, not from the previous to the next.
fr_N_ab=95;
Stats_crt=TRACKING.(sprintf('frame_%d', fr_N_ab));
Stats_crt(Stats_crt.Centroid_p(:,1)==0 | Stats_crt.Centroid_x(:,1)==0,:)=[];
figure
imshow(importedvideo(:,:,fr_N_ab-seq_adj)*2+importedvideo(:,:,fr_N_ab-seq_adj+1),[0,5])
hold on
quiver(Stats_crt.Centroid(:,1),Stats_crt.Centroid(:,2),Stats_crt.Centroid_x(:,1)-Stats_crt.Centroid(:,1),Stats_crt.Centroid_x(:,2)-Stats_crt.Centroid(:,2),'LineWidth',2,'color','yellow','AutoScale','off')
axis on
axis ij
hold off
str_title=['Trajectory from Frame ',num2str(fr_N_ab),' to Frame ',num2str(fr_N_ab+1)];
title(str_title)