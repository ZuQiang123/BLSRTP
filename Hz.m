function  y = Hz(X,Y,Z,x,y,z,I0,D0)

L = cos(I0)*cos(D0);
M = cos(I0)*sin(D0);
N = sin(I0);


R =  sqrt((X-x).^2+(Y-y).^2+(Z-z).^2);

Bx = 2*atan(((y-Y)+R+(x-X))./(z-Z));

By = log((x-X)+R);

Bz = log((y-Y)+R);


y = (Bx*N + By*M + Bz*L);

end
% MMM Magnetization inversion method


