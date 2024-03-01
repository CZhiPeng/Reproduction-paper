%% HF-SVD-ZF,全连接
function [R_HF] = HF_ZF_SVD_fullyconnected(Nt,Num_users,H,H1,Wrf,wk,SNR)
          %速率初始化
          R_HF=0;

         for u=1:Num_users
            %模拟预编码Frf
            Hk = H( :,:,u);% user channel of i  %czp:H = Nt*Nr*K
            tk=wk(:,u)'*Hk;
            tk_angle=angle(tk); %1*Nt
            Frf(:,u)=(1./sqrt(Nt))*exp(1i*(-1)*tk_angle);
         end


        %HF-SVD的有效模拟通道(effective analog channel)
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
            R_HF=R_HF+log2(det(1+SINRk_HF));
        end
end

