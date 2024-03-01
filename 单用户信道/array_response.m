%输入方位角和俯仰角、天线个数
%输出导向矢量。注：已经归一化

% a1是方位角，a2是俯仰角，N是天线总数(sqrt(N)*sqrt(N))
function y = array_response(a1,a2,N)
for m= 0:sqrt(N)-1
    for n= 0:sqrt(N)-1
        y(m*(sqrt(N))+n+1) = exp( 1i* pi* ( m*sin(a1)*sin(a2) + n*cos(a2) ) );
    end
end
y = y.'/sqrt(N);  %y.'表示y的转置，y'表示y的共轭转置
end