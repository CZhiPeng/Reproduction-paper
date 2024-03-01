function [R_BF] = BF_ZF(Num_users,H,H1,Wrf,wk,SNR)
        %% BF-ZF
        R_BF=0;
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
            R_BF=R_BF+log2(det(1+SINRk_BF));
        end
end

