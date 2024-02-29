function  [UserVar,rho,rhow,g]=DefineDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B)

fprintf('setting Bedmachine densities \n');
    
rho = 917;
rhow = 1027;

g = 9.81/1000;

fprintf('Done loading densities \n');
