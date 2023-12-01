function T = prism_m(x,y,z,xi1,xi2,yi1,yi2,zi1,zi2,I,D,I1,D1)
%% ¼ÆËã¡÷T×Ü³¡

temp1 = myfun2(x,y,z,xi2,yi2,zi2,I,D,I1,D1);
temp2_1= myfun2(x,y,z,xi1,yi2,zi2,I,D,I1,D1);
temp2_2 = myfun2(x,y,z,xi2,yi1,zi2,I,D,I1,D1);
temp2_3= myfun2(x,y,z,xi2,yi2,zi1,I,D,I1,D1);
temp3_1 = myfun2(x,y,z,xi2,yi1,zi1,I,D,I1,D1);
temp3_2 = myfun2(x,y,z,xi1,yi2,zi1,I,D,I1,D1);
temp3_3 = myfun2(x,y,z,xi1,yi1,zi2,I,D,I1,D1);
temp4 = myfun2(x,y,z,xi1,yi1,zi1,I,D,I1,D1);
T = (temp1-temp2_1-temp2_2-temp2_3+temp3_1+temp3_2+temp3_3-temp4)*100;