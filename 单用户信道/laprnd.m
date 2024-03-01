function y  = laprnd(m, n, mu, sigma)
%LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
%生成独立同分布从拉普拉斯分布中提取的拉普拉斯随机数
%   with mean mu and standard deviation sigma. 
%   mu      : mean  %平均值
%   sigma   : standard deviation  %标准差
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1. 
%   For more information, refer to
%   http://en.wikipedia.org./wiki/Laplace_distribution

%   Author  : Elvis Chen (bee33@sjtu.edu.cn)
%   Date    : 01/19/07

%Check inputs
if nargin < 2
    error('At least two inputs are required'); %这里检查输入参数的数量。如果提供的输入参数数量小于2，将引发错误，因为至少需要提供 m 和 n 这两个参数。
end

if nargin == 2
    mu = 0; sigma = 1; %如果输入参数的数量等于2，将使用默认值 mu = 0 和 sigma = 1，这表示生成标准拉普拉斯分布的随机数。
end

if nargin == 3 %如果输入参数的数量等于3，将使用默认值 sigma = 1，这表示生成具有自定义均值 mu 的拉普拉斯分布的随机数。
    sigma = 1;
end

% Generate Laplacian noise  %翻译：生成拉普拉斯信号
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
% y = mu - b.* sign(u).* log(1- 2.* abs(u));
y = mu + sign(u) .* ( pi - b.* log( exp(pi./b) + (2-2.*exp(pi./b)) .* abs(u) ) );