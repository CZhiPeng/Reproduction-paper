function [Wrf,wk] = Gain_Wrf(Nr,Num_users,H)
%GAIN_WRF 此处显示有关此函数的摘要
%   此处显示详细说明
        for u=1:Num_users

            %用户k的信道
            Hk = H( :,:,u);% user channel of i  %czp:H = Nt*Nr*K

            %用户k的结合矩阵wk
            [U,~,~]=svd(Hk);
            u1=U(:,1);
            u1_angle=angle(u1); %Nr*1
            wk(:,u)=(1./sqrt(Nr))*exp(1i*u1_angle);
            wk_cell{u} = wk(:, u);
        end

        Wrf=blkdiag(wk_cell{:}); 

end

