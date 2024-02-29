function [trans]=SSTREAM_Tus_t_3d_m_symb(kx,ky,t,alph,H,eta,C,rho,g)

syms m
assume(m,'real');

% time-dependent ratio between surface topography and
% bedrock in Fourier space in dimensional units

% for non-dimensional a la gudmundsson 2003
% put eta=1/2, H0=1,rho*g=1/sin(alpha), lambda in units of mean ice thickness
% Cnondimensional=2 eta/H0 Cdimesional, and then
% U0nondimensional=Cnondimensional 

[k,l] = kxky2kl(kx,ky);


j2=k.^2+l.^2;
ca=cot(alph);
tau=rho*g*sin(alph)*H;
%U=C*tau^m;

gamm=tau^(1-m) / (C*m); % can't use gamma because it is a function
ps=m.*gamm.^2+H.*(l.^2.*(4+m)+k.^2.*(1+4.*m)).*gamm.*eta+4.*H.^2.*j2.^2.*eta.^2;

%upsilon=ca.*H.*k.*gamm;
%xi=gamm+4*H*j2*eta;


%t1=expp.*tau.*(gamm+1i.*ca.*H^2.*k.*j2*eta+H.*(k.^2+4.*l.^2).*eta+1i.*upsilon);
%t2=H.*(gamm+H.*j2.*eta).*xi;
 
t1=(m.*(gamm+1i.*ca.*H.*k.*gamm)+H.*(k.^2.*(1+1i.*ca.*H.*k)+(4+1i.*ca.*H.*k).*l.^2).*eta).*tau;
t2=H.*ps;


trans = abs(t1./t2);

trans = simplify(abs(trans))

return

end



