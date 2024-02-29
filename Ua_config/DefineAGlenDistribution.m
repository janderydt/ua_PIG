function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
    
persistent FAGlen

if ~exist(CtrlVar.NameOfFileForReadingAGlenEstimate)
    AGlen=s*0+AGlenVersusTemp(-15);
    n=3;
    fprintf('\n Using constant intial rate factor %s \n',num2str(AGlen(1)));
else
    if isempty(FAGlen)
        load(CtrlVar.NameOfFileForReadingAGlenEstimate,'xA','yA','AGlen');
        FAGlen = scatteredInterpolant(xA,yA,AGlen,'linear');
        fprintf('\n Read rate factor from file %s \n',CtrlVar.NameOfFileForReadingAGlenEstimate);
    end
    
    load(CtrlVar.NameOfFileForReadingAGlenEstimate,'n');
    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
    AGlen=FAGlen(x,y);
    n = n(1);
end
end
