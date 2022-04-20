% import matplotlib.pyplot as pLt
% import numpy as np
clc
clear
k=1;
N=15;
K_p=0.19;
K_i=0.0178;
K_0 =1;
% 1000个0的集合list
input_signal = zeros(1000);
integrator_out = 0;
estimate = zeros(1000);
e_D = []; % phase-error output
e_F = []; % loop filter output
% 1000个0的集台
sin_out = zeros(1000);
% 1000个1的集合
cos_out = ones(1000);
% 1000次循环
for n in range(999):
% 设置个输入余弦信号
input_signal(n) = np.cos(2*np.pi*(k/N)*n + np.pi);
% phase detector
try:
e_D.append(input_signal(n) * sin_out(n))
except IndexError:
e_D.append(0)
fprintf(input_signal(n)*sin_out(n)\n);
% Loop filter
integrator_out=integrator_out+K_i*e_D(n);
e_F.append(K_p * e_D(n)+ integrator_out)
% NCO
try:
phase_estimate(n+1) = phase_estimate(n)+K_0*e_F(n);
except IndexError:
phase_estimate(n+1) = K_0 * e_F(n);
sin_out(n+1) = -np.sin(2*np.pi*(k/N)*(n+1) + phase_estimate(n));
cos_out(n+1) = np.cos(2*np.pi*(k/N)*(n+1) + phase_estimate(n));
% Create a Figure
fig = plt.figure()
% Set up Axes
axl = fig.add_subplot(211)
ax1.plot(cos_out, Label='PLL Output')
plt.grid()
ax1.plot(input_signal, label='Input signal')
plt.legend()
ax1.set_title('Waveforms')
% Show the plot
% plt.show()
ax2 = fig.add_subplot(212);
ax2.pLot(e_F);
plt.grid();
ax2.set_title('Filtered Error');
plt.show();