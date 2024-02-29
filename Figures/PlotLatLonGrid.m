function [Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing,Colour,lineWidth, LabelFontSize)

%
% Plots a lat lon grid
%[Lat,Lon,X0,Y0,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing,Colour,LineWidth,LabelFontSize)
%
%           scale :  distance units, if using m set to 1 (default), if using km set to 1000
%      dlat, dlon :  spacing between lat and lon lines (default is dlat=2, dlon=5)
%          Colour : color of lines and labels, e.g. [0.8 0.8 0.8] for light grey
%       lineWidth : line width of the lat/lon lines
%   LabelFontSize : font size of the labels
%
% on return
%
% Lat, Lon the meshgrid with lat lon values
% Clat and Clon are the contour matrixes
% hlat and hlon the contour objects
%
% Example:
% [~,~,~,~,Clat,hlat,Clon,hlon]=PlotLatLonGrid(scale,dlat,dlon,LabelSpacing)
% 
% To plot a Lat and Lon-grid for the entire Antarctic continent:
%   PlotLatLonGrid(1000, 5, 30, 2, [0.8 0.8 0.8], 0.1, 4);
%

tt=axis;
xmin=tt(1) ; xmax=tt(2) ; ymin=tt(3) ; ymax=tt(4) ;

if nargin<1 || isempty(scale)
    scale=1;
end

if nargin<2 || isempty(dlat)
    dlat=2;
end

if nargin<3 || isempty(dlon)
    dlon=5;
end

if nargin<4 || isempty(LabelSpacing)
    LabelSpacing =2;
end


if nargin<5 || isempty(Colour)
   Colour='black'; 
end

if nargin<6 || isempty(lineWidth)
    lineWidth = 0.1;
end

if nargin<7 || isempty(LabelFontSize)
    LabelFontSize = 4;
end

lcol='k';

% create lat/lon variables
[X0,Y0]=meshgrid(linspace(xmin,xmax,400),linspace(ymin,ymax,400));
[Lat,Lon]=pol_to_geog_wgs84_71S(X0*scale,Y0*scale);


hold on

% plot latitude lines
levels = [-90:dlat:0];
[Clat,hlat]=contour(X0,Y0,Lat,levels,'--','LineColor',lcol,'LineWidth',lineWidth);
vlat = levels(2:LabelSpacing:end); % set manually, if the labels are not convenient
clabel(Clat,hlat,vlat,'Color','k','FontSize',LabelFontSize,'LabelSpacing',400);

% do not plot longitudes for lat <-85 
Lon(Lat<-85)=NaN;

Lon2 = Lon; % for line at -180=180

% avoid many height lines between 180 and -180 by setting the region between -175 and -180 to NaN 
I=X0<0 & Y0<0 & Lon<-175;
Lon(I)=NaN;

% plot longitudes
levels=-180+dlon:dlon:180;
[Clon,hlon]=contour(X0,Y0,Lon,levels,'--','LineColor',lcol,'LineWidth',lineWidth);
vlon = levels(1:LabelSpacing:end); % set manually, if not convenient, e.g. [ -90  0 90 180];
clabel(Clon,hlon,vlon,'Color','k','FontSize',LabelFontSize,'LabelSpacing',400);
%%% add a line at -180=180 degree %%%

% set everything bigger that -170 and smaller than 170 to NaN
I = Lon2 >-170 & Lon2 < 170;
Lon2(I) = NaN;

% set everythin below -170 to 190 
I = Lon2 < -170;
Lon2(I) = 190;

[Clon2,hlon2]=contour(X0,Y0,Lon2,levels,'LineColor',lcol,'LineWidth',lineWidth);
vlon2 = [180];
clabel(Clon2,hlon2,vlon2,'FontSize',LabelFontSize);


%%% set colors of labels and lines %%%

hlon.LineColor  = Colour ; 
hlon2.LineColor = Colour ;
hlat.LineColor  = Colour ; 
clabel(Clat,  hlat, 'Color',Colour)
clabel(Clon,  hlon, 'Color',Colour)
clabel(Clon2, hlon2,'Color',Colour)

end
