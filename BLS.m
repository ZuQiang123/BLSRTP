%%
%建立训练集及网络训练的程序  
%在本研究中我们假设没有剩余磁化强度
%训练样本正演公式参照袁洋等发表在地球物理学报上的《基于 BTTB 矩阵的快速高精度三维磁场正演》
%我们把研究区域划分成若干网格，计算每个网格在观测点上产生的磁异常
%输入表示倾斜磁场
%标签表示垂直磁场
%网络选用宽度学习网络，参考陈俊龙等发表的《Broad Learning System: An Effective and Efficient Incremental
% Learning System Without the Need for Deep Architecture》

%Code for creating test data
%We assume no remanent magnetization in this study
%The forward  formula for creating training samples refers to "Fast and high accuracy 3D magnetic anomaly
% forward modeling based on BTTB matrix" published by Yuan Yang et al. in Chinese Journal of Geophysics (in Chinese)
%We divide the study area into several grids and calculate the magnetic anomalies (TMI) produced by each grid at the observation points.
%Input means the oblique TMI
%Label means the vertical TMI
%Network selects the breadth learning network, refer to 《Broad Learning System: An Effective and Efficient Incremental
%Learning System Without the Need for Deep Architecture》published by Chen Junlong et al.


clear
clc
tic

infile = 'Test_Input.grd';                                    %读取测试数据  Load test Input
[Ta,x,y,z,nx,ny,dx,dy] = readgrd(infile);         %读取测试数据  Load test Input
%%
% 建立正演矩阵   building sensitivity matrix

Xc = x;  Yc = y;                                                   % 网格中心坐标   Center coordinates of grids
Zc = 0.5:9.5;   dz = 1;                                         % 网格中心坐标   Center coordinates of grids in Z dicrecion   
xi1 = Xc-dx/2; xi2 = Xc+dx/2;                           %网格起始位置   X direction grid boundary
yi1 = Yc-dy/2; yi2 = Yc+dy/2;                           %网格起始位置   Y direction grid boundary
zi1 = Zc-dz/2; zi2 = Zc+dz/2;                           %网格起始位置   Z direction grid boundary

a = length(Xc);    b = length(Yc);   c = length(Zc);
l = length(Xc);     v = length(Yc);

Xtext = reshape(Ta,1,[]);                                 %倾斜磁场，即训练样本的输出     Input, oblique TMI                                

mNumber = b*a*c;    nNumber = l*v;
G = zeros(l*2,v*2,c);                                        %输入核函数    sensitivity matrix
Gh = zeros(l*2,v*2,c);                                      %标签核函数    sensitivity matrix

I = deg2rad(0);                                                 %输入的磁化倾角  Inclination for Input+
D = deg2rad(0);                                               %输入的磁化偏角  Declination for Input
Ih = deg2rad(90);                                             %标签的磁化倾角  Inclination for Label
Dh = deg2rad(0);                                              %标签的磁化偏角  Inclination for Label

for k=1:c                                                           %正演矩阵的计算，可查看参考文献
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
Num = 20000;                                                    %训练样本数量  Number of training samples
Xtrain = zeros(Num,a*b);                                   %训练样本输入 Input of Training samples
Ytrain = zeros(Num,a*b);                                   %训练样本输出 Label of Training samples
for  iNum=1:Num

    m = zeros(2*a,2*b,c);                                     %网格的磁化强度大小  Magnetization(A/m) for each grids

    x_bin1 = randi([1,a],1,1);                                %正演模型的随机坐标  Random coordinates of the forward model
    x_bin2 = randi([x_bin1,a],1,1);                       %正演模型的随机坐标  Random coordinates of the forward model
    y_bin1 = randi([1,b],1,1);                               %正演模型的随机坐标  Random coordinates of the forward model
    y_bin2 = randi([y_bin1,b],1,1);                       %正演模型的随机坐标  Random coordinates of the forward model
    z_bin1 = randi([1,c/2],1,1);                            %正演模型的随机坐标  Random coordinates of the forward model
    z_bin2 = randi([z_bin1,c],1,1);                       %正演模型的随机坐标  Random coordinates of the forward model
    MagnNum = randi([1,50],1,1);                       %磁化强度取值范围     Magnetization intensity value range

    for k=z_bin1:z_bin2                                       %随机建立单个长方体正演模型
        for j = y_bin1:y_bin2                                  %Randomly build a single cuboid forward model
            for i =x_bin1:x_bin2
                m(i,j,k) =  MagnNum;
            end
        end
    end


    TaT =zeros(a*2,b*2);
    Thh =zeros(a*2,b*2);

    for i=1:c

        TaT = ifft2(fft2(G(:,:,i)).*fft2(m(:,:,i)))+TaT;
        Thh = ifft2(fft2(Gh(:,:,i)).*fft2(m(:,:,i)))+Thh;
    end


    Ta = TaT(1:a,1:b);                                %倾斜磁场，即输出     Input, oblique TMI
    Th = Thh(1:a,1:b);                                %垂直磁场，即标签     Label, oblique TMI
    Xtrain(iNum,:) = reshape(Ta,1,[]);
    Ytrain(iNum,:) = reshape(Th,1,[]);

end

%% BLS Training

train_x = (Xtrain);
train_y = (Ytrain);
test_x = (Xtext);
addpath (genpath('BLS_regression')) % add relative path
addpath (genpath('Dispersion_curves')) % add relative path
addpath (genpath('Normalization')) % add relative path
train_y_add = train_y;
for i = 1:size(train_y)
    temp = train_y(i,:);
    index = find(temp == min(temp));
    train_y_add(i,index) = (rand(1,1)*ones(1,length(index))-1)*20;
end
train_x_add = zeros(size(train_x,1),size(train_x,2))+9999;
train_x = [train_x;train_x_add];
train_y = [train_y;train_y_add];
[mean_x,std_x,train_x_norm] = normalized_fun(train_x);
[mean_y,std_y,train_y_norm] = normalized_fun(train_y);
test_x_norm = zeros(size(test_x,1),size(test_x,2));
for i = 1:1:size(test_x,2)
    test_x_norm(:,i) = (test_x(:,i)-mean_x(i))/std_x(i);
end


%%

Fea_vec = 200;                                          %特征节点神经元数目       Number of feature node neurons
Win_vec = 100;                                          %特征节点数目                 Number of feature nodes
Enhan_vec = 100;                                      %增强节点数目                 Number of Enhancement nodes

[Y_hat_norm,NumFea_hat,NumWin_hat,NumEnhan_hat] = ...
    bls_regression_Y_sub_withTrain(train_x_norm,train_y_norm,Fea_vec,Win_vec,Enhan_vec,test_x_norm);
Y_hat = zeros(size(Y_hat_norm,1),size(Y_hat_norm,2));

for i = 1:1:size(train_y_norm,2)
    Y_hat(:,i) = (Y_hat_norm(:,i)*std_y(i)) + mean_y(i);
end

toc
%%


Ttrain = zeros(a,b);

aaa = 0;
for j=1:b
    for i=1:a
        aaa = aaa+1;

        Ttrain(i,j) = Y_hat(1,aaa);

    end
end

figure()
contourf(Xc,Yc,Ttrain)%%画图查看训练数据
colorbar
grid on
xlabel('W-E (m)')
ylabel('S-N (m)')
colormap("jet")










