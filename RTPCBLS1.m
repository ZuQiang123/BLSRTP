clear
clc
tic



bian = 1;
%% 
% 建立正演矩阵
x = 1:30; y=1:30; z = -0.001;

dx = 1;  dy = 1;  dz = 1;% 立方体边长
Xc = x;  Yc = y;   Zc = 0.5:11.5;%  立方体中心坐标
I = 60;    D = 0;      %磁化强度方向



xi1 = Xc-dx/2; xi2 = Xc+dx/2;%积分起始位置
yi1 = Yc-dy/2; yi2 = Yc+dy/2;%积分起始位置
zi1 = Zc-dz/2; zi2 = Zc+dz/2;%积分起始位置

a = length(Xc);    b = length(Yc);   c = length(Zc);
l = length(Xc);     v = length(Yc);


mNumber = b*a*c;    nNumber = l*v;
G = zeros(l*2,v*2,c);%核函数

for k=1:c

    for i=1:l
        for j=1:v

            G(i,j,k) = prism_m(x(i),y(j),z,xi1(1),xi2(1),...
                yi1(1),yi2(1),zi1(k),zi2(k),I,D,I,D);

        end
    end

    for i=1:l
        for j=1:v-1

            G(i,j+v+1,k) = prism_m(x(i),y(1),z,xi1(1),xi2(1),...
                yi1(v-j+1),yi2(v-j+1),zi1(k),zi2(k),I,D,I,D);

        end
    end

    for i=1:l-1
        for j=1:v

            G(i+l+1,j,k) = prism_m(x(1),y(j),z,xi1(l-i+1),xi2(l-i+1),...
                yi1(1),yi2(1),zi1(k),zi2(k),I,D,I,D);
        end
    end

    for i=1:l-1
        for j=1:v-1

            G(i+l+1,j+v+1,k) = prism_m(x(1),y(1),z,xi1(l-i+1),xi2(l-i+1),...
                yi1(v-j+1),yi2(v-j+1),zi1(k),zi2(k),I,D,I,D);
        end
    end

end

%%
load('Xtext1.mat')
iNum=0;
NN = 0;
Num = 10000;
Xtrain = zeros(Num,a*b);
Ytrain = zeros(Num,a*c*b);
while  iNum<Num
    NN = NN+1;
    m = zeros(2*a,2*b,c);

    x_bin1 = randi([10,15],1,1);
    x_bin2 = randi([16,20],1,1);
    y_bin1 = randi([7,15],1,1);
    y_bin2 = randi([16,30],1,1);
    z_bin1 = randi([1,c/2],1,1);
    z_bin2 = randi([z_bin1,c],1,1);
    MagnNum = randi([1,50],1,1);   %磁化强度取值范围

    for k=z_bin1:z_bin2
        for j = y_bin1:y_bin2
            for i =x_bin1:x_bin2
                m(i,j,k) =  MagnNum/10;
            end
        end
    end

    TaT =zeros(a*2,b*2);


    for i=1:c

        TaT = ifft2(fft2(G(:,:,i)).*fft2(m(:,:,i)))+TaT;

    end

    Ta = TaT(1:a,1:b);
    Xtextb(1,:) = reshape(Ta,1,[]);

    AA = corrcoef(Xtext1,Xtextb);

    if AA(2,1)>0.9
        iNum = iNum+1;
        Xtrain(iNum,:) = Xtextb;
        Ytrain(iNum,:) = reshape(m(1:a,1:b,:),1,[]);
    end

end
toc
%%



Ttrain = zeros(a,b);
aaa = 0;
for j=1:b
    for i=1:a
        aaa = aaa+1;

        Ttrain(i,j) = Xtrain(650,aaa);

    end
end
AA = corrcoef(Xtext1,Xtrain(650,:));
aaaa = string(AA(2,1));
figure()
contourf(Ttrain)%%画图查看训练数据
colorbar
text(24,28,aaaa,"FontSize",14,"Color",'r')


load('Ytext1.mat')
Ttext = zeros(a,b);
aaa = 0;
for j=1:b
    for i=1:a
        aaa = aaa+1;

         Ttext (i,j) =  Xtext1(1,aaa);

    end
end
figure()
contourf(Ttext)%%画图查看训练数据
colorbar


GLG = zeros(1,a*b);
aaa = 0;
for j=1:b
    for i=1:a
        aaa = aaa+1;
         GLG(1,aaa) =Ttext(i,j);

    end
end


train_x = (Xtrain);
train_y = (Ytrain);
test_x = (GLG);


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

Fea_vec = 200;
Win_vec = 100;
Enhan_vec = 100;

[Y_hat_norm,NumFea_hat,NumWin_hat,NumEnhan_hat] = ...
    bls_regression_Y_sub_withTrain(train_x_norm,train_y_norm,Fea_vec,Win_vec,Enhan_vec,test_x_norm);
Y_hat = zeros(size(Y_hat_norm,1),size(Y_hat_norm,2));
for i = 1:1:size(train_y_norm,2)
    Y_hat(:,i) = (Y_hat_norm(:,i)*std_y(i)) + mean_y(i);
end

Y_hat(Y_hat<0)=0;

toc

aaa = 0;
M = zeros(a,b,c);
for k = 1:c
    for j=1:b
        for i=1:a
            aaa = aaa+1;
            M(i,j,k) = Y_hat(aaa);
        end
    end
end

aaa = 0;
Mtext = zeros(a,b,c);
for i = 1:a
    for j=1:b
        for k=1:c
            aaa = aaa+1;
            Mtext(i,j,k) = Ytext1(aaa);
        end
    end
end



m = zeros(c,b);
mt = zeros(c,b);
for i=1:c
    for j=1:b

        m(i,j) = M(15,j,i);
        mt(i,j) = Mtext(15,j,i);
    end
end

figure()
subplot(2,1,1)
imagesc(m)
daspect([1 1 1])
colorbar

grid on
subplot(2,1,2)
imagesc(mt)
daspect([1 1 1])
colorbar

grid on

M1 = M;
save('M1','M1')


aaa = 0;
m = zeros(a*b*c,1);
mt = zeros(a*b*c,1);
for i=1:a
    for j=1:b
        for k=1:c
            aaa = aaa+1;
            m(aaa) = M(i,j,k);
            mt(aaa) = Mtext(i,j,k);
        end
    end
end


AA = corrcoef(m,mt);

R = AA(2,1)*AA(2,1);

disp(R)














