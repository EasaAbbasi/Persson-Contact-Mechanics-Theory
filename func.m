function [volt,A,ubar,dg]= func(sigma,sigma0,Es)

epsil0 = 8.854187817e-12;
h0 = .3e-6;

u = logspace(-100,100,5000);                % second input should be infinity
[P_usg,A,~] = perss2(u,sigma,sigma0,Es);

ii = []; I = 0;
for i = 1:length(u)
    ii = [ii, P_usg(i)/(u(i)+h0)^2];
    %     jj = [jj, u(i)*P_usg(i)];
end
I = trapz(u,ii);

%ubar = trapz(u,u.*P_usg);
ubar = 1/sqrt(I)-h0;  % more acurate than previous line

% f = sigma - sigma0 - epsil0/2*volt^2*I;
volt = sqrt(2*(sigma - sigma0)/(epsil0*I));

z = logspace(-9,0,100); jj = []; I =0;
for j = 1:length(z)
    ii=[];
    for i = 1:length(u)
        ii = [ii, P_usg(i)/(u(i)+h0+z(j))^2];
        %     jj = [jj, u(i)*P_usg(i)];
    end
    I = trapz(u,ii);
    jj = [jj,I];
end
J = trapz(z,jj);
dg = epsil0*volt^2/2*J;



