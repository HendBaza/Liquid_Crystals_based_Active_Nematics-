% 1 micrometer = 2.5 pixels
% 1 pixel = 0.4 micrometer
% using mass center from imageJ
load('data_2023_01_03_imageJ.mat')
%%
tosave=0;
weighting=1;
fr_rel=0:7; %averaging
for fr_st=61%[11 19 61 111 161 211 261 301 351 401 451 501 551]
    IS=ImagJStats.(sprintf('frame_%d', fr_st));
    data_table=IS(:,[5,6,9,27]);
    cut_r_edge=8;%cut inward of the original frame
    cut_c_edge=12;
    data_table(:,1:2)=data_table(:,1:2)*2.5+[1,1];
    data_table(data_table(:,1)<=cut_c_edge,:)=[];
    data_table(data_table(:,2)<=cut_r_edge,:)=[];
    data_table(data_table(:,1)>1824-cut_c_edge,:)=[];
    data_table(data_table(:,2)>1216-cut_r_edge,:)=[];
    angs=data_table(:,4); % positive angles are clockwise from x-axis; y-axis points downwards
    angs(angs>pi/2)=pi-angs(angs>pi/2); %adjust angle range from -pi/2 to pi/2
    angs=-angs-pi/2;% positive angles are counterclockwise from y-axis, pointing upwards.
    angs(angs<-pi/2)=pi+angs(angs<-pi/2);%adjust angle range from -pi/2 to pi/2
    data_table(:,4)=angs;
    %%
    for t=16%width of the ring
        t2=t/2;
        kmin=t2+1;
        s1=1216-2*cut_r_edge-1;
        s2=1824-2*cut_c_edge-1;
        kmax=round(min(s1,s2)/4,0);
        for fr_n=fr_rel
            fr_ab=fr_n+fr_st;
            disp('Calculating Frame:')
            disp(fr_ab)
            tic
            Ns=zeros(kmax,1);
            Ns_p=Ns; %for partial rings
            Ws=Ns; %weighted number of bacteria in the ring
            Ws_p=Ns;
            Wms=Ns;
            As=Ns;
            As_ctr=As;
            As_rel2=As_ctr;
            Asoo=As_ctr;
            Asoo_ctr=As_ctr;
            % s_x
            Stats=data_table(data_table(:,3)==fr_n,:);
            Stats(:,1:2)=Stats(:,1:2)-[cut_c_edge,cut_r_edge];% adjust the center region such that (8,12), the left corner is shifted to (0,0)
            Orientation=Stats(:,4);
            Op=sum(exp(2i*Orientation));
            N=size(Stats,1);
            op=norm(Op)/N;
            opx=sum(cos(2*Orientation))/N;
            rho=N/((s1)*(s2));
            C0=4*pi*rho*t2;
            C2=4/3*pi*rho*t;
            Distances=zeros(N);%squared distances with rows for mass in the center and columns for mass in the whole cut-off frame
            Ang=Distances;%angles pointing from the center of ring to mass center in the ring.
            Ang_ctr=Ang;
            Ang_rel2=Ang;
            for i=1:N-1
                xi=Stats(i,1);
                yi=Stats(i,2);
                for j=i+1:N
                    xj=Stats(j,1);
                    yj=Stats(j,2);%y-coordinates are up side down.
                    v=[xj-xi,yj-yi];
                    d=norm(v);
                    Distances(i,j)=d;
                    Distances(j,i)=d;
                    ang=atan2(v(2),v(1));%positive angles are clockwise from the x-axis.
                    Ang(i,j)=-ang-pi/2;%positive angles are counterclockwise from the y-axis.(shear direction)
                    Ang(j,i)=-ang+pi/2;
                    ang_rel=-ang-Stats(i,4);%angle between v(line 66) and bacterium orientation
                    ang_rel2=Stats(j,4)-Stats(i,4);
                    if ang_rel>pi/2
                        ang_rel=ang_rel-pi;
                    end
                    if ang_rel<-pi/2
                        ang_rel=ang_rel+pi;
                    end
                    if ang_rel2>pi/2
                        ang_rel2=ang_rel2-pi;
                    end
                    if ang_rel2<-pi/2
                        ang_rel2=ang_rel2+pi;
                    end
                    Ang_ctr(i,j)=ang_rel;
                    Ang_ctr(j,i)=ang_rel+pi;
                    Ang_rel2(i,j)=ang_rel2;
                    Ang_rel2(j,i)=ang_rel2+pi;
                end
            end
            %
            K = (kmin:kmax)';
            if weighting==0
                Wms(kmin:end)=C0*K;
            else
                Wms(kmin:end)=C2*K;
            end

            check=0;
            for i =1:N %i-th bacterium
                if check==1
                    break
                end
                xi=Stats(i,1);
                yi=Stats(i,2);
                ki=min([xi-t2,s2-t2-xi,yi-t2,s1-t2-yi,kmax]);%max ring radiu for i-th bacterium

                for k = kmin:ki %whole rings
                    Ns(k)=Ns(k)+1;
                end

                for k = kmin:kmax %partial rings
                    r_in=k-t2;
                    r_out=k+t2;
                    caseid=-1;
                    part_check=0;
                    if xi>=r_out && s2-xi>=r_out && yi>=r_out && s1-yi>=r_out %No cutting, Case 0
                        Ns_p(k)=Ns_p(k)+1;
                        caseid=0;
                        part_check=part_check+1;
                    else
                        area_ring=pi*(r_out^2-r_in^2);
                        if xi>=r_out && xi+r_out<=s2 %Case 1s
                            if yi<=r_in %Case 1a
                                yib=yi;
                                caseid=1.1;
                                part_check=part_check+1;
                            end
                            if s1-yi<=r_in %Case 1b
                                yib=s1-yi;
                                caseid=1.2;
                                part_check=part_check+1;
                            end
                            if yi>r_in && yi<r_out %Case 1c
                                yib=yi;
                                caseid=1.3;
                                part_check=part_check+1;
                            end
                            if s1-yi>r_in && s1-yi<r_out %Case 1d
                                yib=s1-yi;
                                caseid=1.4;
                                part_check=part_check+1;
                            end
                            temp1=sqrt(r_out^2-yib^2);
                            temp3=acos(yib/r_out);
                            if caseid==1.1 || caseid==1.2
                                temp2=sqrt(r_in^2-yib^2);
                                temp4=acos(yib/r_in);
                                A=0.5*yib*(temp1-temp2)-0.5*r_in^2*(temp3-temp4); %A1
                            end
                            if caseid==1.3 || caseid==1.4
                                A=0.5*yib*temp1-0.5*r_in^2*temp3; %A4
                            end
                            Ns_p(k)=Ns_p(k)+1-acos(yib/r_out)/pi+2*A/area_ring;
                        end
                        if yi>=r_out && s1-yi>=r_out
                            if xi<=r_in %Case 2a
                                caseid=2.1;
                                xib=xi;
                                part_check=part_check+1;
                            end
                            if s2-xi<=r_in %Case 2b
                                caseid=2.2;
                                xib=s2-xi;
                                part_check=part_check+1;
                            end
                            if xi>r_in && xi<r_out %Case 2c
                                xib=xi;
                                caseid=2.3;
                                part_check=part_check+1;
                            end
                            if s2-xi>r_in && s2-xi<r_out %Case 2d
                                xib=s2-xi;
                                caseid=2.4;
                                part_check=part_check+1;
                            end
                            temp1=sqrt(r_out^2-xib^2);
                            temp3=acos(xib/r_out);
                            if caseid==2.1 || caseid==2.2
                                temp2=sqrt(r_in^2-xib^2);
                                temp4=acos(xib/r_in);
                                A=0.5*xib*(temp1-temp2)-0.5*r_in^2*(temp3-temp4); %A2
                            end
                            if caseid==2.3 || caseid==2.4
                                A=0.5*xib*temp1-0.5*r_in^2*temp3; %A2b
                            end
                            Ns_p(k)=Ns_p(k)+1-acos(xib/r_out)/pi+2*A/area_ring;
                        end
                        if xi^2+yi^2<=r_in^2%Case 3s xi<=r_in && yi<=r_in &&
                            caseid=3.1;
                            xib=xi;
                            yib=yi;
                            part_check=part_check+1;
                        end
                        if (s2-xi)^2+yi^2<=r_in^2 %s2-xi<=r_in && yi<=r_in &&
                            caseid=3.2;
                            xib=s2-xi;
                            yib=yi;
                            part_check=part_check+1;
                        end
                        if (s2-xi)^2+(s1-yi)^2<=r_in^2 %s2-xi<=r_in && s1-yi<=r_in &&
                            caseid=3.3;
                            xib=s2-xi;
                            yib=s1-yi;
                            part_check=part_check+1;
                        end
                        if xi^2+(s1-yi)^2<=r_in^2 %xi<=r_in && s1-yi<=r_in &&
                            caseid=3.4;
                            xib=xi;
                            yib=s1-yi;
                            part_check=part_check+1;
                        end
                        if caseid>=3 && caseid<4
                            temp1=sqrt(r_out^2-yib^2);
                            temp3=acos(yib/r_out);
                            temp2=sqrt(r_in^2-yib^2);
                            temp4=acos(yib/r_in);
                            A1=0.5*yib*(temp1-temp2)-0.5*r_in^2*(temp3-temp4); %A1
                            temp11=sqrt(r_out^2-xib^2);
                            temp31=acos(xib/r_out);
                            temp21=sqrt(r_in^2-xib^2);
                            temp41=acos(xib/r_in);
                            A2=0.5*xib*(temp11-temp21)-0.5*r_in^2*(temp31-temp41); %A2
                            Ns_p(k)=Ns_p(k)+(A1+A2)/area_ring+0.5*(1.5-(temp3+temp31)/pi);
                        end
                        if caseid==-1
                            if xi^2+yi^2>r_in^2 && xi^2+yi^2<=r_out^2 && xi<=r_in && yi>r_in && yi<r_out %Case 4a
                                caseid=4.1;
                                xib=xi;
                                yib=yi;
                                part_check=part_check+1;
                            end
                            if (s2-xi)^2+yi^2>r_in^2 && (s2-xi)^2+yi^2<=r_out^2 && s2-xi<=r_in && yi>r_in && yi<r_out %Case 4b
                                caseid=4.2;
                                xib=s2-xi;
                                yib=yi;
                                part_check=part_check+1;
                            end
                            if (s2-xi)^2+(s1-yi)^2>r_in^2 && (s2-xi)^2+(s1-yi)^2<=r_out^2 && s2-xi<=r_in && s1-yi>r_in && s1-yi<r_out %Case 4c
                                caseid=4.3;
                                xib=s2-xi;
                                yib=s1-yi;
                                part_check=part_check+1;
                            end
                            if xi^2+(s1-yi)^2>r_in^2 && xi^2+(s1-yi)^2<=r_out^2 && xi<=r_in && s1-yi>r_in && s1-yi<r_out %Case 4d
                                caseid=4.4;
                                xib=xi;
                                yib=s1-yi;
                                part_check=part_check+1;
                            end
                            if caseid>=4 && caseid<5
                                temp3=acos(yib/r_out);
                                temp11=sqrt(r_out^2-xib^2);
                                temp31=acos(xib/r_out);
                                temp21=sqrt(r_in^2-xib^2);
                                temp41=acos(xib/r_in);
                                A2=0.5*xib*(temp11-temp21)-0.5*r_in^2*(temp31-temp41);
                                A3=(2*yib-sqrt(r_in^2-xib^2))*xib/2-r_in^2/2*asin(xib/r_in);
                                A4=yib*sqrt(r_out^2-yib^2)/2-r_in^2/2*temp3;
                                Ns_p(k)=Ns_p(k)+(A2+A3+A4)/area_ring+0.5*(1.5-(temp3+temp31)/pi);
                            end
                        end
                        if xi^2+yi^2>r_in^2 && xi^2+yi^2<r_out^2 && xi>r_in && xi<r_out && yi<=r_in %Case 5s
                            caseid=5.1;
                            xib=xi;
                            yib=yi;
                            part_check=part_check+1;
                        end
                        if (s2-xi)^2+yi^2>r_in^2 && (s2-xi)^2+yi^2<r_out^2 && s2-xi>r_in && s2-xi<r_out && yi<=r_in
                            caseid=5.2;
                            xib=s2-xi;
                            yib=yi;
                            part_check=part_check+1;
                        end
                        if (s2-xi)^2+(s1-yi)^2>r_in^2 && (s2-xi)^2+(s1-yi)^2<r_out^2 && s2-xi>r_in && s2-xi<r_out && s1-yi<=r_in
                            caseid=5.3;
                            xib=s2-xi;
                            yib=s1-yi;
                            part_check=part_check+1;
                        end
                        if xi^2+(s1-yi)^2>r_in^2 && xi^2+(s1-yi)^2<r_out^2 && xi>r_in && xi<r_out && s1-yi<=r_in
                            caseid=5.4;
                            xib=xi;
                            yib=s1-yi;
                            part_check=part_check+1;
                        end
                        if caseid>=5 && caseid<6
                            temp1=sqrt(r_out^2-yib^2);
                            temp3=acos(yib/r_out);
                            temp2=sqrt(r_in^2-yib^2);
                            temp4=acos(yib/r_in);
                            temp31=acos(xib/r_out);
                            A1=0.5*yib*(temp1-temp2)-0.5*r_in^2*(temp3-temp4);
                            A5=yib*(2*xib-sqrt(r_in^2-yib^2))/2-r_in^2/2*asin(yib/r_in);
                            A6=0.5*xib*sqrt(r_out^2-xib^2)-r_in^2/2*temp31;
                            Ns_p(k)=Ns_p(k)+(A1+A5+A6)/area_ring+0.5*(1.5-(temp3+temp31)/pi);
                        end
                        if xi<=r_in && yi<=r_in && xi^2+yi^2>=r_out^2%Case 6s
                            xib=xi;
                            yib=yi;
                            caseid=6.1;
                            part_check=part_check+1;
                        end
                        if s2-xi<=r_in && yi<=r_in && (s2-xi)^2+yi^2>=r_out^2
                            xib=s2-xi;
                            yib=yi;
                            caseid=6.2;
                            part_check=part_check+1;
                        end
                        if s2-xi<=r_in && s1-yi<=r_in && (s2-xi)^2+(s1-yi)^2>=r_out^2
                            xib=s2-xi;
                            yib=s1-yi;
                            caseid=6.3;
                            part_check=part_check+1;
                        end
                        if xi<=r_in && s1-yi<=r_in && xi^2+(s1-yi)^2>=r_out^2
                            xib=xi;
                            yib=s1-yi;
                            caseid=6.4;
                            part_check=part_check+1;
                        end
                        if caseid>=6 && caseid<7
                            temp1=sqrt(r_out^2-yib^2);
                            temp3=acos(yib/r_out);
                            temp2=sqrt(r_in^2-yib^2);
                            temp4=acos(yib/r_in);
                            temp11=sqrt(r_out^2-xib^2);
                            temp31=acos(xib/r_out);
                            temp21=sqrt(r_in^2-xib^2);
                            temp41=acos(xib/r_in);
                            A1=0.5*yib*(temp1-temp2)-0.5*r_in^2*(temp3-temp4);
                            A2=0.5*xib*(temp11-temp21)-0.5*r_in^2*(temp31-temp41);
                            Ns_p(k)=Ns_p(k)+(2*A1+2*A2)/area_ring+(1-(temp3+temp31)/pi);
                        end
                        if caseid==-1
                            if xi<=r_in && yi<=r_in && xi^2+yi^2>=r_in^2 && xi^2+yi^2<r_out^2 %Case 7s
                                xib=xi;
                                yib=yi;
                                caseid=7.1;
                                part_check=part_check+1;
                            end
                            if s2-xi<=r_in && yi<=r_in && (s2-xi)^2+yi^2>=r_in^2 && (s2-xi)^2+yi^2<r_out^2
                                xib=s2-xi;
                                yib=yi;
                                caseid=7.2;
                                part_check=part_check+1;
                            end
                            if s2-xi<=r_in && s1-yi<=r_in && (s2-xi)^2+(s1-yi)^2>=r_in^2 && (s2-xi)^2+(s1-yi)^2<r_out^2
                                xib=s2-xi;
                                yib=s1-yi;
                                caseid=7.3;
                                part_check=part_check+1;
                            end
                            if xi<=r_in && s1-yi<=r_in && xi^2+(s1-yi)^2>=r_in^2 && xi^2+(s1-yi)^2<r_out^2
                                xib=xi;
                                yib=s1-yi;
                                caseid=7.4;
                                part_check=part_check+1;
                            end
                            if caseid>=7 && caseid<8
                                temp1=sqrt(r_out^2-yib^2);
                                temp3=acos(yib/r_out);
                                temp2=sqrt(r_in^2-yib^2);
                                temp4=acos(yib/r_in);
                                temp11=sqrt(r_out^2-xib^2);
                                temp31=acos(xib/r_out);
                                temp21=sqrt(r_in^2-xib^2);
                                temp41=acos(xib/r_in);
                                A1=0.5*yib*(temp1-temp2)-0.5*r_in^2*(temp3-temp4);
                                A2=0.5*xib*(temp11-temp21)-0.5*r_in^2*(temp31-temp41);
                                A9=(yib-sqrt(r_in^2-xib^2))*xib/2;
                                A10=(xib-sqrt(r_in^2-yib^2))*yib/2;
                                Ns_p(k)=Ns_p(k)+(A1+A2+A9+A10)/area_ring+0.5*(1.5-(temp3+temp31)/pi);
                            end
                        end
                        if caseid==-1
                            if xi>=r_in && yi>=r_in && xi<=r_out && yi<=r_out %Case 8s
                                caseid=8.1;
                                xib=xi;
                                yib=yi;
                                part_check=part_check+1;
                            end
                            if s2-xi>=r_in && yi>=r_in && s2-xi<=r_out && yi<=r_out %Case 8s
                                caseid=8.2;
                                xib=s2-xi;
                                yib=yi;
                                part_check=part_check+1;
                            end
                            if s2-xi>=r_in && s1-yi>=r_in && s2-xi<=r_out && s1-yi<=r_out %Case 8s
                                caseid=8.3;
                                xib=s2-xi;
                                yib=s1-yi;
                                part_check=part_check+1;
                            end
                            if xi>=r_in && s1-yi>=r_in && xi<=r_out && s1-yi<=r_out %Case 8s
                                caseid=8.4;
                                xib=xi;
                                yib=s1-yi;
                                part_check=part_check+1;
                            end
                            if caseid>=8.1 && caseid<9
                                temp1=sqrt(r_out^2-yib^2);
                                temp3=acos(yib/r_out);
                                temp11=sqrt(r_out^2-xib^2);
                                temp31=acos(xib/r_out);
                                B1=r_out^2*temp3-yib*temp1;
                                B2=r_out^2*temp31-xib*temp11;
                                Ns_p(k)=Ns_p(k)+1-(B1+B2)/area_ring;
                            end
                        end
                        if caseid==-1 %Case 9s and 10s
                            if xi^2+yi^2>=r_out
                                if xi<r_in && r_in<yi && yi<r_out %double cut on x %Case 9a
                                    caseid=9.1;
                                    part_check=part_check+1;
                                    xib=xi;
                                    yib=yi;
                                end
                                if yi<r_in && r_in<xi && xi<r_out %double cut on y %Case 9a
                                    caseid=10.1;
                                    part_check=part_check+1;
                                    xib=xi;
                                    yib=yi;
                                end
                            end
                            if (s2-xi)^2+yi^2>=r_out
                                if s2-xi<r_in && r_in<yi && yi<r_out
                                    caseid=9.2;
                                    part_check=part_check+1;
                                    xib=s2-xi;
                                    yib=yi;
                                end
                                if yi<r_in && r_in<s2-xi && s2-xi<r_out
                                    caseid=10.2;
                                    part_check=part_check+1;
                                    xib=s2-xi;
                                    yib=yi;
                                end
                            end
                            if (s2-xi)^2+(s1-yi)^2>=r_out
                                if s2-xi<r_in && r_in<s1-yi && s1-yi<r_out
                                    caseid=9.3;
                                    part_check=part_check+1;
                                    xib=s2-xi;
                                    yib=s1-yi;
                                end
                                if s1-yi<r_in && r_in<s2-xi && s2-xi<r_out
                                    caseid=10.3;
                                    part_check=part_check+1;
                                    xib=s2-xi;
                                    yib=s1-yi;
                                end
                            end
                            if xi^2+(s1-yi)^2>=r_out
                                if xi<r_in && r_in<s1-yi && s1-yi<r_out
                                    caseid=9.4;
                                    part_check=part_check+1;
                                    xib=xi;
                                    yib=s1-yi;
                                end
                                if s1-yi<r_in && r_in<xi && xi<r_out
                                    caseid=10.4;
                                    part_check=part_check+1;
                                    xib=xi;
                                    yib=s1-yi;
                                end
                            end
                            temp1=sqrt(r_out^2-yib^2);
                            temp3=acos(yib/r_out);
                            temp11=sqrt(r_out^2-xib^2);
                            temp31=acos(xib/r_out);
                            if caseid>=9.1 && caseid<10
                                temp21=sqrt(r_in^2-xib^2);
                                temp41=acos(xib/r_in);
                                A2=0.5*xib*(temp11-temp21)-0.5*r_in^2*(temp31-temp41);
                                B1=r_out^2*temp3-yib*temp1;
                                Ns_p(k)=Ns_p(k)+1-temp31+(2*A2-B1)/area_ring;
                            end
                            if caseid>=10.1 && caseid<11
                                temp2=sqrt(r_in^2-yib^2);
                                temp4=acos(yib/r_in);
                                A1=0.5*yib*(temp1-temp2)-0.5*r_in^2*(temp3-temp4);
                                B2=r_out^2*temp31-xib*temp11;
                                Ns_p(k)=Ns_p(k)+1-temp3+(2*A1-B2)/area_ring;
                            end
                        end
                        if caseid==-1
                            if xi^2+yi^2>r_in^2 && xi^2+yi^2<r_out^2 && xi>r_in && xi<r_out && yi>r_in && yi<r_out %Case 11s
                                caseid=11.1;
                                xib=xi;
                                yib=yi;
                                part_check=part_check+1;
                            end
                            if (s2-xi)^2+yi^2>r_in^2 && (s2-xi)^2+yi^2<r_out^2 && s2-xi>r_in && s2-xi<r_out && yi>r_in && yi<r_out
                                caseid=11.2;
                                xib=s2-xi;
                                yib=yi;
                                part_check=part_check+1;
                            end
                            if (s2-xi)^2+(s1-yi)^2>r_in^2 && (s2-xi)^2+(s1-yi)^2<r_out^2 && s2-xi>r_in && s2-xi<r_out && s1-yi>r_in && s1-yi<r_out
                                caseid=11.3;
                                xib=s2-xi;
                                yib=s1-yi;
                                part_check=part_check+1;
                            end
                            if xi^2+(s1-yi)^2>r_in^2 && xi^2+(s1-yi)^2<r_out^2 && xi>r_in && xi<r_out && s1-yi>r_in && s1-yi<r_out
                                caseid=11.4;
                                xib=xi;
                                yib=s1-yi;
                                part_check=part_check+1;
                            end
                            if caseid>11 && caseid<12
                                B3=xib*yib-pi*r_in^2/4;
                                temp1=sqrt(r_out^2-yib^2);
                                temp3=acos(yib/r_out);
                                temp2=sqrt(r_in^2-yib^2);
                                temp4=acos(yib/r_in);
                                temp31=acos(xib/r_out);
                                A4=yib*sqrt(r_out^2-yib^2)/2-r_in^2/2*temp3;
                                A6=0.5*xib*sqrt(r_out^2-xib^2)-r_in^2/2*temp31;
                                Ns_p(k)=Ns_p(k)+(B3+A4+A6)/area_ring+0.5*(1.5-(temp3+temp31)/pi);
                            end
                        end
                    end
                    if part_check~=1 || imag(Ns_p(k))~=0
                        disp('Check')
                        disp(i)
                        disp(k)
                        disp(caseid)
                        disp(check)
                        check=1;
                        break
                    end
                end
                if imag(Ns_p(k))~=0
                    disp('warning')
                    break
                end
                for j=1:N%j-th bacterium
                    if j~=i
                        dij=Distances(i,j);
                        for k2=kmin:kmax %partial
                            d=abs(k2-dij);
                            if dij<=k2+t2 && dij>=k2-t2%d<=t2-0.1
                                if weighting==0
                                    wk=1;
                                else
                                    wk=1-4*d^2/t^2;
                                end
                                Ws_p(k2)=Ws_p(k2)+wk;
                            end
                        end
                        for k2=kmin:ki %whole
                            d=abs(k2-dij);
                            if dij<=k2+t2 && dij>=k2-t2%d<=t2-0.1
                                phij=Ang(i,j);
                                phij_ctr=Ang_ctr(i,j);
                                phij_rel2=Ang_rel2(i,j);
                                if weighting==0
                                    wk=1;
                                else
                                    wk=1-4*d^2/t^2;
                                end
                                Ws(k2)=Ws(k2)+wk;
                                As(k2)=As(k2)+wk*exp(2i*phij);%total anistrophic number of bacteria
                                As_ctr(k2)=As_ctr(k2)+wk*exp(2i*phij_ctr);
                                As_rel2(k2)=As_rel2(k2)+wk*exp(2i*phij_rel2);
                                Asoo(k2)=Asoo(k2)+wk*real(exp(2i*phij_rel2))*exp(2i*phij);
                                Asoo_ctr(k2)=Asoo_ctr(k2)+wk*real(exp(2i*phij_rel2))*exp(2i*phij_ctr);
                            end
                        end
                    end
                end
            end
            %
            WNS_p=Wms.*Ns_p;
            P_p=Ws_p./WNS_p;
            Cdd0_p=P_p-1;
            WNS=Wms.*Ns;
            P=Ws./WNS;
            Cdd0=P-1;
            Cdd2=As./WNS;
            Cdd2_p=As./WNS_p;
            Cdd2_ctr=As_ctr./WNS; %(discarded) the anistrophic part using orientation of the center bacterium
            Coo0=As_rel2./WNS-op^2;
            Coo0x=As_rel2./WNS-opx^2;
            Coo2=Asoo./WNS;
            Coo2_ctr=Asoo_ctr./WNS;
            Ak_frames.(sprintf('frame_%d', fr_ab))=As;
            Ak_rel_frames.(sprintf('frame_%d', fr_ab))=As_ctr;
            Ak_rel2_frames.(sprintf('frame_%d', fr_ab))=As_rel2;
            Asoos_frames.(sprintf('frame_%d', fr_ab))=Asoo;
            Asoos_ctr_frames.(sprintf('frame_%d', fr_ab))=Asoo_ctr;
            Wk_frames.(sprintf('frame_%d', fr_ab))=Ws;
            Wk_frames_p.(sprintf('frame_%d', fr_ab))=Ws_p;
            Wk_m_frames.(sprintf('frame_%d', fr_ab))=Wms;
            Nk_frames.(sprintf('frame_%d', fr_ab))=Ns;
            Nk_frames_p.(sprintf('frame_%d', fr_ab))=Ns_p;
            Cdd0s.(sprintf('frame_%d', fr_ab))=Cdd0;
            Cdd0s_p.(sprintf('frame_%d', fr_ab))=Cdd0_p;
            Cdd2s.(sprintf('frame_%d', fr_ab))=Cdd2;
            Cdd2s_p.(sprintf('frame_%d', fr_ab))=Cdd2_p;
            Cdd2s_ctr.(sprintf('frame_%d', fr_ab))=Cdd2_ctr;
            Coo0s.(sprintf('frame_%d', fr_ab))=Coo0;
            Coo0sx.(sprintf('frame_%d', fr_ab))=Coo0x;
            Coo2s.(sprintf('frame_%d', fr_ab))=Coo2;
            Coo2s_ctr.(sprintf('frame_%d', fr_ab))=Coo2_ctr;
            ops.(sprintf('frame_%d', fr_ab))=op;
            opsx.(sprintf('frame_%d', fr_ab))=opx;
            Ops.(sprintf('frame_%d', fr_ab))=Op;
            N_bact.(sprintf('frame_%d', fr_ab))=N;
            toc
        end
        if tosave==1
            save('data_2023_01_06_Corr_5C0_t16_imgj.mat','N_bact','ops','Ops','Ak_frames','Ak_rel_frames','Ak_rel2_frames','Wk_frames','Wk_frames_p','Wk_m_frames','Nk_frames','Nk_frames_p','Asoos_frames','Asoos_ctr_frames','Cdd0s','Cdd0s_p','Cdd2s','Cdd2s_ctr','Coo0s','Coo2s','Coo2s_ctr')
        end
    end
end
%% plots
sp1=1;
sp2=1;
figure
cnt=0;
for fr_ab=211%[11 19 61 111 161 211 261 301]%the positive angles are counterclockwise from the y-axis.
    cnt=cnt+1;
    subplot(sp1,sp2,cnt)
    % C=Cdd0s_p.(sprintf('frame_%d', fr_ab));
    C=Cdd2AVE.(sprintf('frame_%d', fr_ab));
    kmax=size(C,1);
    x = linspace(0,kmax*0.4);
    ym= x*0;
    yd = -0.01+x*0;
    yu = .01+x*0;
    hold on
    plot((1:kmax)*0.4,real(C),'LineWidth',2)
%     C=Coo0s.(sprintf('frame_%d', fr_ab));
%     plot((1:kmax)*0.4,real(C),'LineWidth',2)
%     C=Cdd2s.(sprintf('frame_%d', fr_ab));
%     plot((1:kmax)*0.4,real(C),'LineWidth',2)
%     C=Coo2s.(sprintf('frame_%d', fr_ab));
%     plot((1:kmax)*0.4,real(C),'LineWidth',2)
    % C=Cdd2s_ctr.(sprintf('frame_%d', fr_ab));
    % plot((1:kmax)*0.4,real(C),'LineWidth',2)
    % C=Coo2s_ctr.(sprintf('frame_%d', fr_ab));
    % plot((1:kmax)*0.4,real(C),'LineWidth',2)
    plot(x,yu,'r--','DisplayName','y=0.01')
    plot(x,ym,'r--','DisplayName','y=0')
    plot(x,yd,'r--','DisplayName','y=-0.01')
    xlim([0,kmax*0.4])
    xlabel('micrometer','FontSize',20)
    %     ylim([-0.3,0.1])
    hold off
%     legend('C^{(0)}_{dd}','C^{(0)}_{oo}','Location','southeast')
end
%% plots comparison for Cdd0
figure
fr_ab=211;
C=Cdd0s.(sprintf('frame_%d', fr_ab));
kmax=size(C,1);
x = linspace(0,kmax*0.4);
ym= x*0;
yd = -0.01+x*0;
yu = .01+x*0;
hold on
plot((1:kmax)*0.4,real(C),'color','red','LineWidth',2)
C=Cdd0s_p.(sprintf('frame_%d', fr_ab));
plot((1:kmax)*0.4,real(C),'color','magenta','LineWidth',2)
C=Cdd0AVE_p.(sprintf('frame_%d', fr_ab));
plot((1:kmax)*0.4,real(C),'color','blue','LineWidth',2)
plot(x,yu,'r--','DisplayName','y=0.01')
plot(x,ym,'r--','DisplayName','y=0')
plot(x,yd,'r--','DisplayName','y=-0.01')
hold off
xlim([0,kmax*0.4])
xlabel('Distance in $\mu$m','interpreter','latex','FontSize',20)
ylabel('Density-density correlation functions','interpreter','latex','FontSize',20)
%     ylim([-0.3,0.1])
legend('C_1','C_2','C_3')
%% plots comparison for Cdd2,Coo0,Coo2
figure
hold on
for fr_ab=[19,61,161,211,351,451,551]%the positive angles are counterclockwise from the y-axis.
    % C=Cdd0AVE.(sprintf('frame_%d', fr_ab));
    C=Cdd2AVE.(sprintf('frame_%d', fr_ab));
    % C=Coo0AVE.(sprintf('frame_%d', fr_ab));
    % C=Coo2AVE.(sprintf('frame_%d', fr_ab));
    kmax=size(C,1);
    plot((1:kmax)*0.4,real(C),'LineWidth',2)
end
 x = linspace(0,kmax*0.4);
    ym= x*0;
    yd = -0.01+x*0;
    yu = .01+x*0;
    plot(x,yu,'r--','DisplayName','y=0.01')
    plot(x,ym,'r--','DisplayName','y=0')
    plot(x,yd,'r--','DisplayName','y=-0.01')
    hold off
    xlabel('Distance ($\mu m$)','interpreter','latex','FontSize',20)
    % ylabel('Average isotropic orientation-orientation correlation function','interpreter','latex','FontSize',20)
    ylabel('Average anisotropic density-density correlation function','interpreter','latex','FontSize',20)
    xlim([0,kmax*0.4])
    %     ylim([-0.3,0.1])
    legend('3 sec','15 sec','28 sec','40 sec','65 sec','88 sec','113 sec','Location','southeast')
%% plots
sp1=2;
sp2=4;
figure
cnt=0;
for fr_ab=551:558%[11 61 111 161 211 261]%the positive angles are counterclockwise from the y-axis.
    cnt=cnt+1;
    subplot(sp1,sp2,cnt)
    C=Cdd0s.(sprintf('frame_%d', fr_ab));
    kmax=size(C,1);
    x = linspace(0,kmax*0.4);
    ym= x*0;
    yd = -0.01+x*0;
    yu = .01+x*0;
    hold on
    plot((1:kmax)*0.4,real(C),'LineWidth',2)
    C=Cdd0s_p.(sprintf('frame_%d', fr_ab));
    plot((1:kmax)*0.4,real(C),'LineWidth',2)
    plot(x,yu,'r--','DisplayName','y=0.01')
    plot(x,ym,'r--','DisplayName','y=0')
    plot(x,yd,'r--','DisplayName','y=-0.01')
    xlim([0,kmax*0.4])
    xlabel('micrometer','FontSize',20)
    %     ylim([-0.3,0.1])
    hold off
    legend('C^{(0)}_{dd}','C^{(0)}_{dd}(partial)','Location','southeast')
end