function northarrow_psxy(X,Y,x0,y0)

[Lat,Lon] = psxy2ll(X,Y,-71,0);
[Lat0,Lon0] = psxy2ll(x0,y0,-71,0);

%find end point of arrow
Lat1 = Lat0 - 1;
Lon1 = Lon0;

[x1,y1] = ll2psxy(Lat1,Lon1,-71,0);

% rescale arrow
Larrow = sqrt((x1-x0).^2+(y1-y0).^2);
m = (y1-y0)/(x1-x0);

XYscale = 1e3;

arrow([x0 y0]/XYscale,[x0+sqrt(1/(m^2+1)) y0+sqrt(m^2/(m^2+1))]/XYscale,15,'color','k','BaseAngle',60,'Ends',2);
text(x0/XYscale,y0/XYscale+5,'N','fontsize',20,'verticalalignment','bottom');



