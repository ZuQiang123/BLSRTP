%%
%建立测试集的程序  
%在本研究中我们假设没有剩余磁化强度
%正演公式参照袁洋等发表在地球物理学报上的《基于 BTTB 矩阵的快速高精度三维磁场正演》
%我们把研究区域划分成若干网格，计算每个网格在观测点上产生的磁异常
%输入表示倾斜磁场
%标签表示垂直磁场

%Code for creating test data
%We assume no remanent magnetization in this study
%The forward  formula refers to "Fast and high accuracy 3D magnetic anomaly
% forward modeling based on BTTB matrix" published by Yuan Yang et al. in Chinese Journal of Geophysics (in Chinese)
%We divide the study area into several grids and calculate the magnetic anomalies (TMI) produced by each grid at the observation points.
%Input means the oblique TMI
%Label means the vertical TMI
%%

clear
clc
tic

% 建立正演矩阵   building sensitivity matrix

x = 1:50; y=1:50; z = -0.001;       %观测点坐标   Observation point coordinates
dx = 1;  dy = 1;  dz = 1;               % 网格边长  Size of grids
Xc = x;  Yc = y;   Zc = 0.5:9.5;     % 网格中心坐标   Center coordinates of grids
xi1 = Xc-dx/2; xi2 = Xc+dx/2;    %网格起始位置   X direction grid boundary
yi1 = Yc-dy/2; yi2 = Yc+dy/2;    %网格起始位置   Y direction grid boundary
zi1 = Zc-dz/2; zi2 = Zc+dz/2;    %网格起始位置   Z direction grid boundary

a = length(Xc);    b = length(Yc);   c = length(Zc);
l = length(Xc);     v = length(Yc);


mNumber = b*a*c;    nNumber = l*v;
G = zeros(l*2,v*2,c);                 %输入核函数    sensitivity matrix
Gh = zeros(l*2,v*2,c);               %标签核函数    sensitivity matrix


I = deg2rad(0);                         %输入的磁化倾角  Inclination for Input
D = deg2rad(0);                       %输入的磁化偏角  Declination for Input
Ih = deg2rad(90);                     %标签的磁化倾角  Inclination for Label
Dh = deg2rad(0);                     %标签的磁化偏角  Inclination for Label

for k=1:c                                  %正演矩阵的计算，可查看参考文献
                                                 %Calculation of sensitivity matrix, please refer to the references
    for i=1:l    
        for j=1:v

            G(i,j,k) = prism_m(x(i),y(j),z,xi1(1),xi2(1),...
                yi1(1),yi2(1),zi1(k),zi2(k),I,D,I,D);
            Gh(i,j,k) = prism_Za(x(i),y(j),z,xi1(1),xi2(1),...
                yi1(1),yi2(1),zi1(k),zi2(k),Ih,Dh);

        end
    end

    for i=1:l
        for j=1:v-1

            G(i,j+v+1,k) = prism_m(x(i),y(1),z,xi1(1),xi2(1),...
                yi1(v-j+1),yi2(v-j+1),zi1(k),zi2(k),I,D,I,D);
            Gh(i,j+v+1,k) = prism_Za(x(i),y(1),z,xi1(1),xi2(1),...
                yi1(v-j+1),yi2(v-j+1),zi1(k),zi2(k),Ih,Dh);


        end
    end

    for i=1:l-1
        for j=1:v

            G(i+l+1,j,k) = prism_m(x(1),y(j),z,xi1(l-i+1),xi2(l-i+1),...
                yi1(1),yi2(1),zi1(k),zi2(k),I,D,I,D);
            Gh(i+l+1,j,k) = prism_Za(x(1),y(j),z,xi1(l-i+1),xi2(l-i+1),...
                yi1(1),yi2(1),zi1(k),zi2(k),Ih,Dh);

        end
    end

    for i=1:l-1
        for j=1:v-1

            G(i+l+1,j+v+1,k) = prism_m(x(1),y(1),z,xi1(l-i+1),xi2(l-i+1),...
                yi1(v-j+1),yi2(v-j+1),zi1(k),zi2(k),I,D,I,D);
            Gh(i+l+1,j+v+1,k) = prism_Za(x(1),y(1),z,xi1(l-i+1),xi2(l-i+1),...
                yi1(v-j+1),yi2(v-j+1),zi1(k),zi2(k),Ih,Dh);
        end
    end

end


%%

m = zeros(2*a,2*b,c);   %网格的磁化强度大小  Magnetization(A/m) for each grids

for k=2:6                      %可由研究者自己设定   Can be set by researchers
    for j = 20:30             %可由研究者自己设定   Can be set by researchers
        for i = 20:30         %可由研究者自己设定   Can be set by researchers
            m(i,j,k) =  1;     %可由研究者自己设定   Can be set by researchers
        end
    end
end


TaT =zeros(a*2,b*2);    %倾斜磁场，即输出     Input, oblique TMI
Thh =zeros(a*2,b*2);    %垂直磁场，即标签     Label, oblique TMI

for i=1:c

    TaT = ifft2(fft2(G(:,:,i)).*fft2(m(:,:,i)))+TaT;  
    Thh = ifft2(fft2(Gh(:,:,i)).*fft2(m(:,:,i)))+Thh;

end

Ta = TaT(1:a,1:b);
Th = Thh(1:a,1:b);
%%

figure()
contourf(Xc,Yc,Ta)
c = colorbar;
title(c,'Magnetic Anomaly data (nT)','rotation',90,'Fontname', 'Times New Roman','FontSize',14,'position',[54.8,128.58,0]) 
grid on
colormap("jet")
set(gca,'YDir','normal');
xlabel('X (m)')
ylabel('Y (m)')
set(gca,'Fontsize',12)

figure()
contourf(Xc,Yc,Th)
c = colorbar;
title(c,'Magnetic Anomaly data (nT)','rotation',90,'Fontname', 'Times New Roman','FontSize',14,'position',[54.8,128.58,0])
grid on
colormap("jet")
set(gca,'YDir','normal');
xlabel('X (m)')
ylabel('Y (m)')
set(gca,'Fontsize',12)

%%
nx = a;
ny = b;
outfile = 'Test_Input.grd'; 
savegrd(Ta,x,y,z,nx,ny,outfile); %保存输入    save Input
outfile = 'Test_Label.grd'; 
savegrd(Th,x,y,z,nx,ny,outfile);  %保存标签    save Label









       










