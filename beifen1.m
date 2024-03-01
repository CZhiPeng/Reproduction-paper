%备份：没有转成函数之前

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


%% main
for i_snr = 1 :nSNR
    i_snr
    SNR =SNR_linear(i_snr);
    for iter=1:iterations


        % Generate user channels 
        [H,At,Ar]=Multi_user_channel_realization(Nt,Nr,Num_users,Nc,Nray,angle_sigma);
        % H is a 3-dimensional matrix, 规模：Nr*Nt*Num_users

        %多用户信道   czp:对三维矩阵H进行处理，将其变为2维矩阵。使其符合论文中(Num_users × Nr) × Nt的多用户信道
        H1=permute(H,[1 3 2]); %调换H的第2列和第3列,此时的维度是Nr*Num_users*Nt
        H1=reshape(H1,[],Nt); %这样处理后，H_eff 的维度是 (Num_users * Nr) * Nt


        %% HF-SVD-ZF,全连接
        %获得 Wrf、Frf
        for u=1:Num_users

            %用户k的信道
            Hk = H( :,:,u);% user channel of i  %czp:H = Nt*Nr*K

            %用户k的结合矩阵wk
            [U,~,~]=svd(Hk);
            u1=U(:,1);
            u1_angle=angle(u1); %Nr*1
            wk(:,u)=(1./sqrt(Nr))*exp(1i*u1_angle);

            %模拟预编码Frf
            tk=wk(:,u)'*Hk;
            tk_angle=angle(tk); %1*Nt
            Frf(:,u)=(1./sqrt(Nt))*exp(1i*(-1)*tk_angle);
        end

        %用户整体的结合矩阵Wrf
        for u=1:Num_users
            wk_cell{u} = wk(:, u);
        end

        %HF-SVD的有效模拟通道(effective analog channel)
        Wrf=blkdiag(wk_cell{:}); 
        Heff_hf=Wrf'*H1*Frf;

        %HF-SVD,基带是ZF预编码
        hf_f_zf=(sqrt(Num_users))/norm((Frf*Heff_hf'*inv(Heff_hf*Heff_hf')),'fro');
        HF_Fbb_ZF=hf_f_zf*Heff_hf'*inv(Heff_hf*Heff_hf');

        %HF-SVD,和速率
        for u_hf=1:Num_users
            %每个用户信号初始化
            Interference_signal_hf=0;

            Hk = H( :,:,u_hf);
            useful_signal_hf=(wk(:,u_hf)'*Hk*Frf*HF_Fbb_ZF(:,u_hf))^2;
            %useful_signal=(norm((wk(:,u)'*Hk*Frf*Fbb_ZF(:,u)),'fro'))^2;

            for uu_hf=1:Num_users  %累加干扰信号
                if uu_hf~=u_hf
                    Interference_signal_hf=Interference_signal_hf+(wk(:,uu_hf)'*Hk*Frf*HF_Fbb_ZF(:,uu_hf))^2;
                    %Interference_signal=Interference_signal+(norm((wk(:,u_u)'*Hk*Frf*Fbb_ZF(:,u_u)),'fro'))^2;
                end
            end

            SINRk_HF=useful_signal_hf/(Interference_signal_hf+(1./(SNR)));
            %R(i_snr)=R(i_snr)+log2(det(1+inv(Interference_signal+1)*useful_signal));
            R_HF(i_snr)=R_HF(i_snr)+log2(det(1+SINRk_HF));
        end

        %% BF-ZF
        %数字预编码的有效模拟通道(effective analog channel)
        Heff_bf=Wrf'*H1;

        %BF-SVD,ZF预编码
        bf_f_zf=(sqrt(Num_users))/sqrt((trace(inv(Heff_bf*Heff_bf'))));
        BF_Fbb_ZF=bf_f_zf*Heff_bf'*inv(Heff_bf*Heff_bf');
        
        %BF 和速率
        for u_bf=1:Num_users
            %每个用户信号初始化
            Interference_signal_bf=0;

            Hk = H( :,:,u_bf);
            useful_signal_bf=(wk(:,u_bf)'*Hk*BF_Fbb_ZF(:,u_bf))^2;

            for uu_bf=1:Num_users  %累加干扰信号
                if uu_bf~=u_bf
                    Interference_signal_bf=Interference_signal_bf+(wk(:,uu_bf)'*Hk*BF_Fbb_ZF(:,uu_bf))^2;
                end
            end

            SINRk_BF=useful_signal_bf/(Interference_signal_bf+(1./(SNR)));
            R_BF(i_snr)=R_BF(i_snr)+log2(det(1+SINRk_BF));
        end
        
        %% HF-SVD-ZF

        


    end
end

R_HF = R_HF/iterations;
R_BF=R_BF/iterations;
LineWidth = 1.5;
MarkerSize = 6;

figure;
plot(SNR_dB, abs(R_HF), 'b-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'R1');
hold on;
plot(SNR_dB, abs(R_BF), 'r-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'R2');
hold off;grid on;

