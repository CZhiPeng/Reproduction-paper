% This code realizes 1000 mmWave channels
clear,clc

Ns = 3; % # of streams

Nc = 5; % # of clusters
Nray = 10; % # of rays in each cluster

Nt = 144; % # of transmit antennas
Nr = 36; % # of receive antennas

angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
%角度扩展为10°，化为弧度，即标准差

gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H,根据 H 的归一化条件

realization = 1000; %迭代次数
count = 0;
for reali = 1:realization
    for c = 1:Nc
        AoD_m = unifrnd(0,2*pi,1,2); %每个簇的平均角度(AOD)，也是laprnd函数中的均值
        AoA_m = unifrnd(0,2*pi,1,2); %unifrnd(a,b) 从具有下部端点 a 和上部端点 b 的连续均匀分布中生成一个随机数。
        
        AoD(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(1),angle_sigma); %AOD的方位角的生成，有Nc*Nray个
        AoD(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(2),angle_sigma); %AOD的俯仰角
        AoA(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(1),angle_sigma);
        AoA(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(2),angle_sigma);
    end
    
    H(:,:,reali) = zeros(Nr,Nt); %初始化信道矩阵，是第reali次迭代的信道矩阵
    for j = 1:Nc*Nray
        At(:,j,reali) = array_response(AoD(1,j),AoD(2,j),Nt); %UPA array response
        Ar(:,j,reali) = array_response(AoA(1,j),AoA(2,j),Nr);
        alpha(j,reali) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1); %第j个子径的复增益系数
        H(:,:,reali) = H(:,:,reali) + alpha(j,reali) * Ar(:,j,reali) * At(:,j,reali)';
    end
    H(:,:,reali) = gamma * H(:,:,reali);
    
    if(rank(H(:,:,reali))>=Ns) %取出最优的模拟预编码、数字预编码矩阵
        count = count + 1;
    
        [U,S,V] = svd(H(:,:,reali));
        Fopt(:,:,reali) = V([1:Nt],[1:Ns]);
        Wopt(:,:,reali) = U([1:Nr],[1:Ns]);
    end
end