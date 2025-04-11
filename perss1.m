function [G,A,ubar] = perss1(zt1,sigma1,sigma0,Es)

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

cq = @(q)min(1e-6.*q.^(-2.*(H+1)),1e-6.*(1e3).^(-2.*(H+1)));
cq3 = @(q)q.^3.*min(1e-6.*q.^(-2.*(H+1)),1e-6.*(1e3).^(-2.*(H+1)));


G_zt = @(zt)pi/4.*(Es).^2.*integral(cq3,qL,zt.*qL);
G = G_zt(zt1);
A_zt_sg = @(sigma,zt)A0.*erf(sigma./(2.*sqrt(G_zt(zt))));
A = A_zt_sg(sigma1,zt1);

% sg_zt_sg = @(sigma,zt) sigma0.*A0./(A_zt_sg(sigma,zt));
% w_qzt = @(q,zt) 1/sqrt(pi.*integral(cq3,zt.*q0,q));
% w_q = @(q) 1/sqrt(pi.*integral(cq3,q0,q));
% P_qzt = @(pp,zt)erf(pp/(2.*sqrt(G_zt(zt))));
% int1 = @(pp,q,zt,sigma)(1./pp).*(gama+3.*(1-gama).*(P_qzt(pp,zt)).^2.*...
%     exp(-pp.*(w_qzt(q,zt))./Es));

sg_zt_sg = @(sigma,zt)sigma0.*A0./(A_zt_sg(sigma1,zt1));

q = logspace(log10(zt1)+logq0,logq1,50);
% pp = linspace(sg_zt_sg(sigma1,zt1),1e6*sg_zt_sg(sigma1,zt1),500);
pp = logspace(log10(sg_zt_sg(sigma1,zt1)),log10(1e6*sg_zt_sg(sigma1,zt1)),15);
I1 = 0; I2 = 0; i1 = []; i2 = [];
for i = 1:length(q)
    w_q = 1/sqrt(pi.*integral(cq3,q0,q(i)));
    w_qzt = 1/sqrt(pi.*integral(cq3,min(zt1.*q0,q(i)),q(i)));
    if (i == 1), w_qzt = Inf; end   %%% this line is used to prevent minus (imaginary part for sqrt) which occures for the i = 1
    i2 = [i2, ((q(i))^2*cq(q(i))*w_q)*I1];
    %     a(i) = max(pi.*integral(cq3,zt1.*q0,q(i)),0);  %%% test line
    i1 = [];
    for j = 1:length(pp)

        %         P_qzt = erf(pp(j)/(2.*sqrt(G_zt(zt1))));

        P_qzt = erf(w_qzt/Es*pp(j));
        int1 = (1./pp(j)).*(gama+3.*(1-gama).*(P_qzt).^2).*exp(-(pp(j).*(w_qzt)./Es)^2);
        i1 = [i1, int1];
        %         a(i,j) = int1;
    end
    I1 = trapz(pp,i1);
    a(i) = ((q(i))^2*cq(q(i))*w_q)*I1;  %%%% test line
end
I2 = trapz(q,i2);
ubar = sqrt(pi)*I2;


end
