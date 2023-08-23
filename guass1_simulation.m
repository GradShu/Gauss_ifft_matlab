clear;close all
format long
%若振幅不对称———>振幅添加噪声
%光谱信号振幅12
A = 24;
A1 = 15;
%中心波长lamda0:770nm；lamda1:850nm
lamda1 = 750e-3;
lamda0 = 770e-3;
%光谱线宽FWHM：15.7nm
delta_lamda = 15.7e-3;
%信号波长lamda
lamda = 0.6:0.001:1;
%生成高斯分布函数
G_lamda = A1*exp(2.77*(-(lamda-lamda1).^2/delta_lamda^2));
F_lamda = A*exp(2.77*(-(lamda-lamda0).^2/delta_lamda^2));
I1_lamda = G_lamda+F_lamda;
% figure(1),subplot(2,2,1);plot(lamda,G_lamda);xlabel('λ(um)','FontSize',12);ylabel('Light Intensity','FontSize',12);
% subplot(2,2,2);plot(lamda,F_lamda);xlabel('λ(um)','FontSize',12);ylabel('Light Intensity','FontSize',12);
% subplot(2,2,3);plot(lamda,I1_lamda);xlabel('λ(um)','FontSize',12);ylabel('Light Intensity','FontSize',12);

%转化为波矢k
sigma = 1./lamda;
%重新排序数据，让横轴按照从小到大的顺序存储
sigma_sorted = sort(sigma);
temp_flamda = zeros(length(sigma),1);
for i = 1:length(sigma)
    TempK = sigma_sorted(i);
    k_index = find(sigma == TempK);
    temp_flamda(i,1) = I1_lamda(k_index);
end
I1_sigma = temp_flamda;
% subplot(2,2,4);plot(sigma_sorted,I1_sigma);
% xlabel('σ(um^{-1})','FontSize',12);ylabel('Light Intensity','FontSize',12);

%相位,设定z0为300um
z0 = 300;
%设置采样点点数为1024？一个周期至少8个采样点
N = 2^16;delta_sigma = 1/(2*N*((770e-3)/16));
sigma =(1:N)*delta_sigma;
% sigma = linspace(sigma_start,sigma_end,N);
%相位
phase = -4*pi*z0*sigma;
%I(sigma) = F(sigma)*exp(j*phase)
sigma0 = 1/lamda0;
sigma1 = 1/lamda1;
% F_sigma = A*exp(2.77*(-(1./sigma-1/sigma0).^2/delta_lamda^2));
G_sigma = A1*exp(2.77*(-(1./sigma-1/sigma1).^2/delta_lamda^2));
F_sigma = A*exp(2.77*(-(1./sigma-1/sigma0).^2/delta_lamda^2));
I1_sigma = G_sigma+F_sigma;
I_sigma = I1_sigma.*exp(1i*phase);
figure(2);subplot(2,1,1);plot(sigma,abs(I_sigma));
xlabel('σ(um^{-1})','FontSize',12);ylabel('Light Intensity','FontSize',12);xlim([1.1 1.5]);
subplot(2,1,2);plot(sigma,phase);
xlabel('σ(um^{-1})','FontSize',12);ylabel('Phase','FontSize',12);xlim([1.1 1.5]);
% subplot(3,1,3);plot(sigma,angle(I_sigma)./pi);
% xlabel('σ(um^{-1})','FontSize',12);ylabel('Phase / \pi','FontSize',12);

%对I_sigma进行傅里叶逆变换
S_z = ifft(I_sigma);
%对逆傅里叶变换后数据处理?将波数K转化为z？
%S_z = fftshift(S_z);
delta_sigma = (sigma(end)-sigma(1))/(N-1);
delta_z = 1/(2*N*delta_sigma);
z = delta_z*linspace(1,N,N);
figure(3);
subplot(2,1,1);plot(z,abs(S_z));
xlabel('z(um))','FontSize',12);ylabel('Light Intensity','FontSize',12);xlim([200 400]);
subplot(2,1,2),plot(z,angle(S_z)./pi);
xlabel('z(um)','FontSize',12);ylabel('Phase','FontSize',12);xlim([298 302]);