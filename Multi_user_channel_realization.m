%  czp:代码参考三部分：
% 1：SIC_based_HP_MU_3D，其中频谱效率公式并非常规，因此只用它的信道部分
% 2.交替优化算法，信道部分是单用户信道
% 3.Low-Complexity Hybrid Precoding in Massive Multiuser MIMO Systems：参考了它的信道部分和SE部分(这里没用到)

% 此函数目的是生成多用户毫米波信道数据，这里的H是3维数组，1：Nr、2：Nt、3：Num_users
%A_t是UPA阵列的导向矢量；A_r是ULA阵列的导向矢量
function [H,At,Ar] = Multi_user_channel_realization(Nt,Nr,Num_users,Nc,Nray,angle_sigma)


% Nc = 1; % 每用户单簇
% Nray = 10; % # of rays in each cluster
% Nt = 144; % # of transmit antennas
% Nr = 36; % # of receive antennas
% Num_users=4; %4个用户
%angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
%角度扩展为10°，化为弧度，即标准差
gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H  %根据 H 的归一化条件


%角度数据,接收端是ULA阵列，因此只有方位角数据；离开角AOD；到达角AOA
for u=1:Num_users
        for c = 1:Nc
            AoD_m = unifrnd(0,2*pi,1,2); %每个簇的平均角度(AOD)，也是laprnd函数中的均值
            AoA_m = unifrnd(0,2*pi,1,2); %unifrnd(a,b) 从具有下部端点 a 和上部端点 b 的连续均匀分布中生成一个随机数。
            AoD(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(1),angle_sigma); %AOD的方位角的生成，有Nc*Nray个
            AoD(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoD_m(2),angle_sigma); %AOD的俯仰角
            AoA(1,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(1),angle_sigma);
            %AoA(2,[(c-1)*Nray+1:Nray*c]) = laprnd(1,Nray,AoA_m(2),angle_sigma);
        end
        H(:,:,u) = zeros(Nr,Nt); 
        for j = 1:Nc*Nray
            At(:,j) = array_response_upa(AoD(1,j),AoD(2,j),Nt); %UPA array response
            %Ar(:,j) = array_response_upa(AoA(1,j),AoA(2,j),Nr); %UPA array response
            Ar(:,j) = array_response_ula(AoA(1,j),Nr);
            alpha(j) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1); %第j个子径的复增益系数
            H(:,:,u) = H(:,:,u) + alpha(j) * Ar(:,j) * At(:,j)';
        end

     H(:,:,u) = gamma * H(:,:,u);

end






