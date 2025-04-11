%% Numerical solutions to Persson's contact mechanics theory in
%Persson, B. N. J. "The dependency of adhesion and friction on electrostatic attraction." The Journal of chemical physics 148, no. 14 (2018).
%
% For the mathematical model, refer to the following paper:
% Easa AliAbbasi, MReza Alipour Sormoli, and Cagatay Basdogan. "Frequency-dependent behavior of electrostatic forces
% between human finger and touch screen under electroadhesion. "IEEE Transactions on Haptics" vol. 15, no. 2, pp. 416-428, 2022.


%*** Attention ***%
% Free use is only permitted by refering to: Easa AliAbbasi, MReza Alipour Sormoli, and Cagatay Basdogan. "Frequency-dependent behavior of electrostatic forces
% between human finger and touch screen under electroadhesion. "IEEE Transactions on Haptics" vol. 15, no. 2, pp. 416-428, 2022.

%%
close all
clear
clc
%% Initializations
sigma0 = 10e3;                      % applied pressure
Es = 10e6;                          % effective elastic modulus
h0 = 0.3e-6;                        % effective thickness
eps0 = 8.854187817e-12;             % permittivity of free space

% ranges
s1 = sigma0+1e3;
s2 = sigma0+10^6.4;
%%
% plot presure, average separation, and A/A0
v1 =(2*(logspace(log10(s1),log10(s2))-sigma0)/eps0*(0+h0)^2).^.5;
figure(1); plot(log10(v1),log10(logspace(log10(s1),log10(s2))-sigma0))
hold on
pause(.3)
v1 =(2*(logspace(log10(s1),log10(s2))-sigma0)/eps0*(.1e-6+h0)^2).^.5;
figure(1); plot(log10(v1),log10(logspace(log10(s1),log10(s2))-sigma0))
hold on
pause(.3)
v1 =(2*(logspace(log10(s1),log10(s2))-sigma0)/eps0*(1000e-9+h0)^2).^.5;
figure(1); plot(log10(v1),log10(logspace(log10(s1),log10(s2))-sigma0))
pause(.3)
hold on

figure(1)
xlabel('log_{10}(voltage) (Volt)')
ylabel('log_{10}(P-P_0) (Pa)')
pos1 = get(gcf,'Position');
lim1 = axis(gca);
hold on

figure(2)
xlabel('log_{10}(voltage) (Volt)')
ylabel('log_{10}(Average Separation) (m)')
pos2 = get(gcf,'Position');
set(gcf,'Position', pos2 - [0,1.23*pos1(4),0,0])
xlim(lim1(1:2))
hold on

figure(3)
xlabel('log_{10}(voltage) (Volt)')
ylabel('log_{10}(A_0/A)')
pos3 = get(gcf,'Position');
set(gcf,'Position', pos3 + [1.08*pos1(3),0,0,0])
hold on


i = 1;v=[]; A = []; ubar = [];
for s = logspace(log10(s1),log10(s2))
    tic
    [v(i),A,ubar,dg]=func(s,sigma0,Es);
    A0 = 1;

    toc
    figure(1); scatter(log10(v(i)),log10(s-sigma0),'filled')
    hold on; grid on;
    figure(2); scatter(log10(v(i)),log10(ubar),'filled')
    hold on; grid on;
    figure(3); scatter(log10(v(i)),log10(A(end-1)/A0),'filled')
    hold on, grid on;
    pause(.3)
    i =i+1;
end
