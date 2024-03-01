
function y = array_response_ula(a1,N) % a1是方位角
for m= 0:N-1
        y(m+1) = exp( 1i* pi* ( m*sin(a1) ) );

end
y = y.'/sqrt(N);  %y.'表示y的转置，y'表示y的共轭转置
end