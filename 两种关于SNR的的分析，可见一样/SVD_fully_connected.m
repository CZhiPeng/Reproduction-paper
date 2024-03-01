clear;clc
%%  ----------------------------- System Parameters -------------------------
Num_users=8; % Number of users
Nt=256; %Number of UPA TX antennas
Nr=16; %Number of ULA(UPA) RX antennas

% ----------------------------- Channel Parameters ------------------------
Nc = 1; % 每用户单簇
Nray = 6; % # of rays in each cluster
angle_sigma = 10/180*pi; %角度扩展为10°，化为弧度，即标准差

% ----------------------------- Simulation Parameters ---------------------
%固定噪声功率，信噪比的一种表达,SNR=P/(Num_users*sigma2)=P/(Num_users*1),P=Num_users*SNR
sigma2=1;
SNR_dB=-25:5:0; SNR_linear=10.^(SNR_dB/10.);   
nSNR = length(SNR_dB);

%固定发射功率为1，信噪比的另一种表达，SNR=P/(Num_users*sigma2)=1/(Num_users*sigma2),sigma2=1./(Num_users*SNR_linear)
% P=1;
% SNR_dB=-25:5:0;SNR_linear=10.^(SNR_dB/10.);   
% nSNR = length(SNR_dB);

iterations=100; % Number of iterations
R1 = zeros(nSNR, 1);
R2 = zeros(nSNR, 1);


%% main
for i_snr = 1 :nSNR
    i_snr
    SNR =SNR_linear(i_snr);
    for iter=1:iterations


        % Generate user channels 
        [H,At,Ar]=Multi_user_channel_realization(Nt,Nr,Num_users,Nc,Nray,angle_sigma);
        % H is a 3-dimensional matrix, 规模：Nr*Nt*Num_users

        %czp:对三维矩阵H进行处理，将其变为2维矩阵。使其符合论文中(Num_users × Nr) × Nt的多用户信道
        H1=permute(H,[1 3 2]); %调换H的第2列和第3列,此时的维度是Nr*Num_users*Nt
        H1=reshape(H1,[],Nt); %这样处理后，H_eff 的维度是 (Num_users * Nr) * Nt


        %% HF-SVD
        %获得 Wrf、Frf
        for u=1:Num_users
            Hk = H( :,:,u);% user channel of i  %czp:H = Nt*Nr*K
            %Wrf
            [U,~,~]=svd(Hk);
            u1=U(:,1);
            u1_angle=angle(u1); %Nr*1
            wk(:,u)=(1./sqrt(Nr))*exp(1i*u1_angle);


            %Frf
            tk=wk(:,u)'*Hk;
            tk_angle=angle(tk); %1*Nt
            Frf(:,u)=(1./sqrt(Nt))*exp(1i*(-1)*tk_angle);
        end
        %Wrf
        for u=1:Num_users
            wk_cell{u} = wk(:, u);
        end
        Wrf=blkdiag(wk_cell{:}); 
        Heff=Wrf'*H1*Frf;

        %获得Fbb,ZF
        f_zf=(sqrt(Num_users))/norm((Frf*Heff'*inv(Heff*Heff')),'fro');
        Fbb_ZF=f_zf*Heff'*inv(Heff*Heff');

        %和速率
        for u=1:Num_users
            %每个用户信号初始化
            %useful_signal=0;
            Interference_signal1=0;
            Interference_signal2=0;

            Hk = H( :,:,u);
            %useful_signal=(SNR/Num_users)*(wk(:,u)'*Hk*Frf*Fbb_ZF(:,u))^2;
            useful_signal1=(SNR)*(wk(:,u)'*Hk*Frf*Fbb_ZF(:,u))^2;
            useful_signal2=(wk(:,u)'*Hk*Frf*Fbb_ZF(:,u))^2;
            %useful_signal=(1/Num_users)*(wk(:,u)'*Hk*Frf*Fbb_ZF(:,u))^2;
            %useful_signal=(norm((wk(:,u)'*Hk*Frf*Fbb_ZF(:,u)),'fro'))^2;

            for u_u=1:Num_users  %累加干扰信号
                if u_u~=u
                    Interference_signal1=Interference_signal1+(SNR)*(wk(:,u_u)'*Hk*Frf*Fbb_ZF(:,u_u))^2;
                    Interference_signal2=Interference_signal2+(wk(:,u_u)'*Hk*Frf*Fbb_ZF(:,u_u))^2;
                    %Interference_signal=Interference_signal+(norm((wk(:,u_u)'*Hk*Frf*Fbb_ZF(:,u_u)),'fro'))^2;
                end
            end

            SINRk1=useful_signal1/(Interference_signal1+sigma2);
            SINRk2=useful_signal2/(Interference_signal2+1./(SNR));
            %R(i_snr)=R(i_snr)+log2(det(1+inv(Interference_signal+1)*useful_signal));
            R1(i_snr)=R1(i_snr)+log2(det(1+SINRk1));
            R2(i_snr)=R2(i_snr)+log2(det(1+SINRk2));
        end

    end
end

R1 = R1/iterations;
R2 = R2/iterations;
LineWidth = 1.5;
MarkerSize = 6;
% figure
% plot(SNR_dB, abs(R), 'k-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)

figure;
plot(SNR_dB, abs(R1), 'b-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'R1');
hold on;
plot(SNR_dB, abs(R2), 'r-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'R2');
hold off;grid on;



