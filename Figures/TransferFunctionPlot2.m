function test2
% for non-dimensional a la gudmundsson 2003
% put eta=1/2, H0=1,rho*g=1/sin(alpha), lambda in units of mean ice thickness
% Cnondimensional=2 eta/H0 Cdimesional, and then
% U0nondimensional=Cnondimensional 
m=[1:50];

g = 9.81*1000;
kx= 10;%:100];
ky = 10;
t = 0;
alpha = 0.01;
H0 = 1;
eta = 1/2;
Cdim = 1e-2;
Cnondim = 2*eta/H0*Cdim;
rho = 1/(sin(alpha)*g);
taud=1;
%rho = 920;

%SSTREAM_Tus_t_3d_m_symb(t,H0,eta,rho,g);
SSTREAM_Tus_t_3d_m_symb(kx,ky,t,alpha,H0,eta,Cnondim,rho,g);
SSTREAM_Tvs_t_3d_m_symb(kx,ky,t,alpha,H0,eta,Cnondim,rho,g);

figure; hold on;

for ii=1:length(m)
    trans_us(ii,:)=abs(SSTREAM_Tus_t_3d_m(kx,ky,t,alpha,H0,eta,Cnondim,rho,g,m(ii)));
    trans_vs(ii,:)=abs(SSTREAM_Tvs_t_3d_m(kx,ky,t,alpha,H0,eta,Cnondim,rho,g,m(ii)));
    %plot(kx,trans(ii,:));
end

psi = 1i*kx*H0*cot(alpha);
j2 = kx^2+ky^2;
g1=Cnondim^(-1)*(1+psi)+eta*H0*(j2*psi+kx^2+4*ky^2);
g3=H0*Cnondim^(-2)+H0^2*Cnondim^(-1)*eta*(4*ky^2+kx^2);
g2=4*H0^3*j2^2*eta^2+eta*H0^2*Cnondim^(-1)*(ky^2+4*kx^2);
transJan = abs(g1.*m./(g2.*m+g3));

m_fit = [1:50];

figure; hold on;
for ii=1:length(kx)
    plot(m,trans_us(:,ii),'o-');
    plot(m,trans_vs(:,ii),'d-');
    [p_quad,S_quad] = polyfit(m,trans_us(:,ii)',2);
    [y_fit_quad,delta_quad] = polyval(p_quad,m_fit,S_quad);
    plot(m_fit,y_fit_quad,'b-')
    plot(m_fit,y_fit_quad+2*delta_quad,'c--',m_fit,y_fit_quad-2*delta_quad,'c--')
    plot(m,transJan,'-+m');
    
    %fexp2 = fit(m',trans_us(:,ii),'exp2');
    options = fitoptions('rat11','lower',[-Inf -Inf -Inf],'upper',[Inf Inf Inf],'startpoint',[1 0 1]);
    [frat11,gof] = fit(m',trans_us(:,ii),'rat11',options);
    frat11
    
    plot(frat11,'k');
end
