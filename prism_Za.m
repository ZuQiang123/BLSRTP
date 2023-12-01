function T = prism_Za(x,y,z,xi1,xi2,yi1,yi2,zi1,zi2,I1,D1)


temp1 = Hz(x,y,z,xi2,yi2,zi2,I1,D1);
temp2_1= Hz(x,y,z,xi1,yi2,zi2,I1,D1);
temp2_2 = Hz(x,y,z,xi2,yi1,zi2,I1,D1);
temp2_3= Hz(x,y,z,xi2,yi2,zi1,I1,D1);
temp3_1 = Hz(x,y,z,xi2,yi1,zi1,I1,D1);
temp3_2 = Hz(x,y,z,xi1,yi2,zi1,I1,D1);
temp3_3 = Hz(x,y,z,xi1,yi1,zi2,I1,D1);
temp4 = Hz(x,y,z,xi1,yi1,zi1,I1,D1);
T = (temp1-temp2_1-temp2_2-temp2_3+temp3_1+temp3_2+temp3_3-temp4)*100;