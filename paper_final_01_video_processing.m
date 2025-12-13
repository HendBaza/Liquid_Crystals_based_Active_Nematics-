% Importing video
N_fr_st=81;%absolute starting frame of the original video
N_fr_ed=450;%absolute ending frame of the original video
n_fr=N_fr_ed-N_fr_st+1;
n_r=1216;
n_c=1824;
video_g=zeros(n_r,n_c,n_fr);
for i = N_fr_st:N_fr_ed
img_g=double(Jan262020_Bsub10to1_DSCGTB_spectrum__13wt_10um_25C_20x_200ms_w1(i).cdata); %<--------------imported video file in RBG
video_g(:,:,i)=img_g(:,:,2);%Choose a channel of RGB
end
%
video_g_f=video_g;
video_g_f(video_g<=80)=NaN;% <80 after shearing
bgmode=ceil(mode(-video_g_f,3));%Putting a negative sign to find max-value mode.
bgstd=round(std(video_g_f,0,3,'Omitnan'),0);
bgstd(isnan(bgstd))=0;
video_g_f(isnan(video_g_f))=0;
%
bgmode(isnan(bgmode))=0;
video_g_nbg=video_g_f+bgmode;
factor=1;%<-----------------------0.9*bgstd before shearing,1*bgstd after shearing
video_g_nbg(video_g_nbg<=factor*bgstd)=0;
% binary videos
video_g_nbg_bi=video_g_nbg;
video_g_nbg_bi(video_g_nbg>0)=1;

video_g_nbg_bi_dn=zeros(size(video_g));
% video_g_nbg_dn=video_g_nbg_bi_dn;
for i=N_fr_st:N_fr_ed
    video_g_nbg_bi_dn(:,:,i)=bwareaopen(video_g_nbg_bi(:,:,i),6);
%     video_g_nbg_dn(:,:,i)=video_g_nbg(:,:,i).*video_g_nbg_bi_dn(:,:,i);
end
%
implay(video_g_nbg_bi_dn)