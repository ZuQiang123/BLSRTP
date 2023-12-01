function  y = myfun2(x,y,z,xi,yi,zi,I,D,I1,D1)
a = length(x);    
k1 = 2*cos((I))*cos((D))...%L
    *cos((I1))*cos((D1));%l
k2 = 2*cos((I))*sin((D))...%M
    *cos((I1))*sin((D1));%m
k3 = 2*sin((I))*sin((I1));%N*n
k4 = sin((I))*cos((I1))*sin((D1))...%N*m
    +cos((I))*sin((D))*sin((I1));%M*n
k5 = sin((I))*cos((I1))*cos((D1))...%N*l
    +cos((I))*cos((D))*sin((I1));%L*n
k6 = cos((I))*sin((D))*cos((I1))*cos((D1))...%M*l
    +cos((I))*cos((D))*cos((I1))*sin((D1));%L*m
temp1 = ones(a,1);      temp4 = ones(a,1);    R = ones(a,1);
temp2 = ones(a,1);      temp5 = ones(a,1);
temp3 = ones(a,1);      temp6 = ones(a,1);
for j = 1:a   
        R(j) = sqrt((x(j)-xi)^2+(y(j)-yi)^2+(z-zi)^2);
        temp1(j) = atan((xi-x(j))/((zi-z)+(yi-y(j))+R(j)));
        temp2(j) = atan((yi-y(j))/((xi-x(j))+(zi-z)+R(j)));
        temp3(j) = atan(((xi-x(j))+(yi-y(j))+R(j))/(zi-z));
        temp4(j) = log((xi-x(j))+R(j));
        temp5(j) = log((yi-y(j))+R(j));
        temp6(j) = log((zi-z)+R(j));  
end
y=-k1.*temp1-k2.*temp2+k3.*temp3+k4.*temp4+k5.*temp5+k6.*temp6;
end



