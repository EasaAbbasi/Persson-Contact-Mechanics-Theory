function [P_usg,AA,uu] = perss2(u,sigma1,sigma0,Es)

E = 100e3;
v = .4;
epsil0 = 8.854187817e-12;
ep_0 = 8.85e-12;
ep_1 = 3.9*ep_0;
ep_2 = 1*ep_0;
ep_3 = 3e3*ep_0;
d_1 = 1e-6;
d_2 = 0.1e-6;
d_3 = 25e-6;
h1 = d_1/ep_1;
h2 = d_3/ep_3;
h0 = h1+h2;
% sigma0 = 50e3; % nominal uniform pressure or stress
H = .71;
hrms = 20e-6;
qL = 1e2;
q0 = 1e3;
q1 = 1e10;
logqL = 2;
logq0 = 3;
logq1 = 10;
A0 = 1;
% Es = E/(1-v^2);
gama = .45;

cqq = @(q)q.*min(1e-6.*q.^(-2.*(H+1)),1e-6.*(1e3).^(-2.*(H+1)));

% zt1 = linspace(1.1,q1/q0,50); dzt = .1;
zt1 = logspace(log10(1.1),log10(floor(q1/q0)),50); dzt = .001;

for i = 1:length(zt1)
    h2rms(i) = integral(cqq,zt1(i)*q0,q1);
    [~,A1(i),ubar1(i)] = perss1(zt1(i) ,sigma1,sigma0,Es);
    [~,A2(i),ubar2(i)] = perss1(zt1(i)+dzt ,sigma1,sigma0,Es);
    A_pr(i) = -(A1(i)-A2(i))/dzt;
    ubar_pr(i) = -(ubar1(i)-ubar2(i))/dzt;
    u1(i) = ubar1(i)+ubar_pr(i)*A1(i)/A_pr(i);
    %     u1(i) = (ubar1(i)*A1(i)-ubar2(i)*A2(i))/(A1(i)-A2(i));
    %     a(i) = h2rms(i);


end
% h2rms(1) = h2rms(2); %%% the same reason for perss1: if (i == 1), w_qzt = 0; end
% a = ubar1/sqrt(h2rms(1));
a = sqrt(h2rms/h2rms(1));
AA = A1;
% uu = ubar1;

kk = 0;
I = []; ii = [];
for j = 1:length(u)
    Ii = 0; ii =[];
    for k = 1:length(zt1)
        %         if (sqrt(2*pi*h2rms(k)) ~= 0)
        kk = kk + 1;
        int2 = (-A_pr(k)/sqrt(2*pi*h2rms(k)))*(exp(-(u(j)-u1(k))^2/(2*h2rms(k)))+exp(-(u(j)+u1(k))^2/(2*h2rms(k))));
        ii = [ii, int2];
        %         a(j,k) = (-A_pr(k)/sqrt(2*pi*h2rms(k)))*...
        %             (exp(-(u(j)-u1(k))^2/(2*h2rms(k)))+exp(-(u(j)+u1(k))^2/(2*h2rms(k))));; %%%
        aa(j,k) = h2rms(k);
        %         end
    end
    %     a(j,:) = ii;
    Ii = trapz(zt1(1:end-1),ii(1:end-1));
    I = [I,Ii];
end


% disp((i-kk))
P_usg = I/A0;
uu = trapz(u,u.*P_usg);
end

