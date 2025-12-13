% load('data_2023_02_07_velocityfield.mat', 'VelocityField_ave','VelocityField_ave_297')
% load('data_2023_01_30_Tracking_5C0_velocityfield.mat', 'TRACKING','TRACKING_296')
%load('data_2023_01_30_tracking_5C0_all.mat')
load('data_velocityfield_2024_01_11_10C0.mat')
%% v=pixel/sec
time=0.4;%sec
savedata=1;
st=90;%11:299%<-------------------------------------------------------------------------------------------------range of frames to be analyzed
nd=101;%304:598;%297:303;%14:296;%
for fr_N_ab=st+1:nd-1%11:299%<-------------------------------------------------------------------------------------------------range of frames
    %fprintf('Analyzing Frame: ')
%     tracking=VelocityField_single.(sprintf('frame_%d', fr_N_ab))(:,1:3);
    Stats_crt_new=TRACKING.(sprintf('frame_%d', fr_N_ab));
    Stats_crt_new(Stats_crt_new.Centroid_p(:,1)==0 | Stats_crt_new.Centroid_x(:,1)==0,:)=[];
    tracking=Stats_crt_new(:,2);
    tracking.DisplaceVel=zeros(size(Stats_crt_new(:,2)));%pixel/sec
    tracking.DisplaceVel=(Stats_crt_new.Centroid_x-Stats_crt_new.Centroid_p)/time;
    tracking.Orientation=Stats_crt_new.Orientation;

    v_mean=mean(tracking.DisplaceVel);
    tracking.DisplaceVel=tracking.DisplaceVel-v_mean;%<----------------------------------------------------mean velocity subtracted
    phi=tracking.Orientation;
    n=[cos(phi),sin(phi)];
    v=tracking.DisplaceVel;
    tracking.Per=v-dot(n,v,2).*n;
    Velocity_no_mean.(sprintf('frame_%d', fr_N_ab))=tracking;
end
%% velocity field
r=48;% Number of the bacteria in the radius of 19.2 microns (48 pixels).
r2=r^2;
for fr_N_ab=st+2:nd-2%15:597%297:303%14:296%11:299%304:596%<-------------------------------------------------------------------------------------------------
    disp(fr_N_ab)
    Stats_crt_new=Velocity_no_mean.(sprintf('frame_%d', fr_N_ab));
    vv=[Stats_crt_new.Orientation,Stats_crt_new.Per];

    Stats_crt_p=Velocity_no_mean.(sprintf('frame_%d', fr_N_ab-1));
    vv_p=[Stats_crt_p.Orientation,Stats_crt_p.Per];

    Stats_crt_x=Velocity_no_mean.(sprintf('frame_%d', fr_N_ab+1));
    vv_x=[Stats_crt_x.Orientation,Stats_crt_x.Per];

    loc=Stats_crt_new.Centroid;
    N=size(Stats_crt_new,1);
    loc_p=Stats_crt_p.Centroid;
    loc_x=Stats_crt_x.Centroid;
    tic
    V_new=table;
    V_new.Centroid=Stats_crt_new.Centroid;
    V_new.DisplaceVel=Stats_crt_new.DisplaceVel;
    V_new.Orientation=Stats_crt_new.Orientation;%ranging from -pi/2 to pi/2. postive angles are from positive x-axis to positive y-axis.
    V_new.Per=Stats_crt_new.Per;
    V_new.Flow=zeros(N,2);
    V_new.Swim=zeros(N,2);
    V_new.BactN=zeros(N,1);% Number of the bacteria in the radius of 19.2 microns (48 pixels).
    for i = 1:N
        loc_s=loc(i,:);
        [M1,V_per1,K1] = tensorform(loc,loc_s,vv,r2);%current
        [M0,V_per0,K0] = tensorform(loc_p,loc_s,vv_p,r2);%previous
        [M2,V_per2,K2] = tensorform(loc_x,loc_s,vv_x,r2);%next
        v_f=[M0;M1;M2]\[V_per0;V_per1;V_per2];
        V_new.Flow(i,:)=v_f;
        V_new.BactN(i)=K0+K1+K2;
    end
    V_new.Swim=V_new.DisplaceVel-V_new.Flow;
    VelocitiesFields_no_mean.(sprintf('frame_%d', fr_N_ab))=V_new;
    toc
end
%%
if savedata==1
    Filename = sprintf('data_%s_10C0_velocitiesfields_no_mean.mat', datestr(now,'yyyy_mm_dd_HH_MM'));
    save(Filename,'VelocitiesFields_no_mean','Velocity_no_mean')
end
%% Draw Fields of velocities with mean velocity removed.
mag1=1;
fr_N_ab=95;%<-----------------------------------------------------------------------------------------------choose a frame to show
V=VelocitiesFields_no_mean.(sprintf('frame_%d', fr_N_ab));
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