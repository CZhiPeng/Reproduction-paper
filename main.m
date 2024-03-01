clear;clc
%%  ----------------------------- System Parameters -------------------------
Num_users=6; % Number of users
Nt=144; %Number of UPA TX antennas
Nr=16; %Number of ULA(UPA) RX antennas

% ----------------------------- Channel Parameters ------------------------
Nc = 1; % 每用户单簇
Nray = 6; % # of rays in each cluster
angle_sigma = 10/180*pi; %角度扩展为10°，化为弧度，即标准差

% ----------------------------- Simulation Parameters ---------------------
%固定噪声功率，信噪比的一种表达,SNR=P/(Num_users*sigma2)=P/(Num_users*1),P=Num_users*SNR
% sigma2=1;
% SNR_dB=-25:5:0; SNR_linear=10.^(SNR_dB/10.);   
% nSNR = length(SNR_dB);

%固定发射功率为1，信噪比的另一种表达，SNR=P/(Num_users*sigma2)=1/(Num_users*sigma2),sigma2=1./(Num_users*SNR_linear)
P=1;
SNR_dB=-30:5:0;SNR_linear=10.^(SNR_dB/10.);   
nSNR = length(SNR_dB);

iterations=100; % Number of iterations
R_BF = zeros(nSNR, 1);
R_HF = zeros(nSNR, 1);
R_HF_Sub= zeros(nSNR, 1);


%% main
for i_snr = 1 :nSNR
    i_snr
    SNR =SNR_linear(i_snr);
    for iter=1:iterations


        % Generate user channels ，H is a 3-dimensional matrix, 规模：Nr*Nt*Num_users
        [H,At,Ar]=Multi_user_channel_realization(Nt,Nr,Num_users,Nc,Nray,angle_sigma);

        %多用户信道，二维   czp:对三维矩阵H进行处理，将其变为2维矩阵。使其符合论文中(Num_users × Nr) × Nt的多用户信道
        H1=permute(H,[1 3 2]); %调换H的第2列和第3列,此时的维度是Nr*Num_users*Nt
        H1=reshape(H1,[],Nt); %这样处理后，H_eff 的维度是 (Num_users * Nr) * Nt

        %结合矩阵，纯模拟
        [Wrf,wk] = Gain_Wrf(Nr,Num_users,H);

        %% HF-SVD-ZF,全连接

        R1=HF_ZF_SVD_fullyconnected(Nt,Num_users,H,H1,Wrf,wk,SNR);
        R_HF(i_snr)=R_HF(i_snr)+R1;

        %% BF-ZF
        R2 = BF_ZF(Num_users,H,H1,Wrf,wk,SNR);
        R_BF(i_snr)=R_BF(i_snr)+R2;

        %% HF-SVD-ZF，半连接
        R3=HF_ZF_SVD_subconnected(Nt,Num_users,H,H1,Wrf,wk,SNR);
        R_HF_Sub(i_snr)=R_HF_Sub(i_snr)+R3;

    end
end

R_HF = R_HF/iterations;
R_BF=R_BF/iterations;
R_HF_Sub=R_HF_Sub/iterations;
LineWidth = 1.5;
MarkerSize = 6;

figure;
plot(SNR_dB, abs(R_BF), 'r-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'BF-ZF');hold on;
plot(SNR_dB, abs(R_HF), 'b-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'HF-SVD-ZF-fullyconnected');
plot(SNR_dB, abs(R_HF_Sub), 'g-d', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'HF-SVD-ZF-subconnected');
hold off;grid on;
legend('Location', 'northwest'); % 设置图例位置在图的左上角
