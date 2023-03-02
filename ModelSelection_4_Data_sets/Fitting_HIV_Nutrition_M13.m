clear all
close all
clc

global tforward initial_cond ViralData CD4Data AlbuminData GlobulinData


ViralData = [8.51E+6	1.13E+7	9.17E+5	3.82E+5	2.73E+5	1.57E+5	1.55E+5	9.51E+4]';
CD4Data = 1000*[1304	619.5	371.1666667	151.6666667	235.8333333	241.1666667	270	379.1666667	246.1666667]';
AlbuminData = [3.44299233423559;2.94258941669938;3.12718763559725;...
               2.92043928342666;3.87677955245160;3.22351592040829;...
               3.51332727235159;3.34054796784925;3.06362997706469]; %mg/dL

GlobulinData = [2.57179227742938;2.59199354586829;2.70640482476377];

data_points = length(ViralData)+length(CD4Data) +...
              length(AlbuminData)+length(GlobulinData) ;


initial_cond = [1304000 0 200 10 3.44 2.57179227742938];

tVLdata = [7,2*7,3*7,4*7,6*7,8*7,10*7,12*7];
tCD4data = [0,7,2*7,3*7,4*7,6*7,8*7,10*7,12*7];
tAlbumindata = [0, 7, 7*2, 7*3,5*7,6*7, 8*7, 12*7, 16*7];  
tGlobulindata = [0, 12*7, 24*7];
tforward = 0:0.1:170;
 
k = [46083.0686392091,0.161003028125379,1.07778443294132e-08,...
     0.896700363149943,2826.05378877325,9.92963416295260,...
     3.30502074055700e-07,8.11958539446260e-08,0.280043973070943,...
     7.67667183288055,1.41187037459318e-08,2.33137194222245,...
     0.999792374260813,3.17704811206496e-09,0.377457207980577,...
     5.26530143560413e-10,0.00277981294517048, 0.1, 0.1];

lb = zeros(1,length(k));

[k,fval] =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter','MaxFunEvals',1e+4, 'MaxIter', 1e+4));


KK = length(k) +1;

AIC = data_points*log(fval/data_points) + 2*KK +(2*KK*(KK+1))/(data_points - KK -1);

[t_r, y_r] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,k),tforward,initial_cond);


display('Parameters after data fitting Model 5:')

fprintf('r = %g\n',  k(1));   
fprintf('d = %g\n', k(2));
fprintf('rho_0 = %g\n', k(3));
fprintf('delta = %g\n',  k(4));
fprintf('p = %g\n', k(5));
fprintf('c = %g\n',  k(6));
fprintf('psi_01 = %g\n', k(7));
fprintf('b = %g\n',  k(8));
fprintf('mu_z = %g\n',  k(9));
fprintf('lambda_A = %g\n', k(10));
fprintf('gamma_A = %g\n', k(11));
fprintf('mu_A = %g\n', k(12));
fprintf('lambda_G = %g\n', k(13));
fprintf('gamma_G = %g\n', k(14));
fprintf('mu_G = %g\n', k(15));
fprintf('mu_v = %g\n', k(16));
fprintf('lambda_z = %g\n', k(17));
fprintf('A_1 = %g\n', k(18));
fprintf('A_2 = %g\n', k(19));
fprintf('AIC = %g\n', AIC);
fprintf('SSE = %g\n', fval);

figure(1)
plot(tforward, y_r(:,3),'-b','LineWidth',2)
hold on 
plot(tVLdata, ViralData, 'ro')
title('Viral Load')

figure(2)
plot(tforward, log10(y_r(:,3)),'-b','LineWidth',2)
hold on 
plot(tVLdata, log10(ViralData), 'ro')
title('Viral Load')

figure(3)

plot(tforward, y_r(:,1),'-b','LineWidth',2)
hold on 
plot(tCD4data, CD4Data, 'ro')
title('CD4 cells')

figure(4)

plot(tforward, log10(y_r(:,1)),'-b','LineWidth',2)
hold on 
plot(tCD4data, log10(CD4Data), 'ro')
title('CD4 cells')

figure(5)

plot(tforward, y_r(:,4),'-b','LineWidth',2)

title('Immune cells')

figure(6)

plot(tforward, y_r(:,5),'-b','LineWidth',2)

title('Albumin')

figure(7)

plot(tforward, y_r(:,6),'-b','LineWidth',2)

title('Globulin')

 function error_in_data = err_in_data(k) 
 
 global  tforward initial_cond ViralData CD4Data AlbuminData GlobulinData %opts = odeset('NonNegative',[1,2,3,4,5,6]);
  
  [~,y] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,k),tforward,initial_cond);
  
  t_v_measure = [7,2*7,3*7,4*7,6*7,8*7,10*7,12*7]./0.1 +1;
  t_cd4_measure = [0,7,2*7,3*7,4*7,6*7,8*7,10*7,12*7]./0.1 +1;
  t_alb_measure = [0, 7, 7*2, 7*3,5*7,6*7, 8*7, 12*7, 16*7]./0.1 +1;
  t_glob_measure = [0,12*7, 24*7]./0.1 +1;
 
 Model_Viral = log10(y(t_v_measure(:),3));
 Model_CD4 = log10(y(t_cd4_measure(:),1));
 Model_Albumin = y(t_alb_measure(:),5); 
 Model_Globulin = y(t_glob_measure(:),6); 
 
 error_in_data = sum((Model_Viral - log10(ViralData)).^2) +...
              20*sum((Model_CD4 - log10(CD4Data)).^2)+...
              sum((Model_Albumin - AlbuminData).^2) +...
              sum((Model_Globulin - GlobulinData).^2); 
 
 end

function dy = Model_HIV_Control_WithinHost(y,k)

dy = zeros(6,1);

r = k(1);
d = k(2);
rho_0 = k(3);
delta = k(4);
p = k(5);
c = k(6);

psi_01 = k(7);
b = k(8);
mu_z = k(9);

lambda_A = k(10); %[15]
gamma_A = k(11);
mu_A = k(12); %[0.33]

lambda_G = k(13); %[0.5,4]
gamma_G = k(14);
mu_G = k(15); %lambda_g/35
mu_v = k(16);
lambda_z = k(17);
A_1 = k(18);
A_2 = k(19);




T = y(1);
T_i = y(2);
V = y(3);
Z = y(4);
A = y(5);
G = y(6);


dy(1) = r - (rho_0/(1 + A_1*A + A_2*G))*T*V - d*T;
dy(2) = (rho_0/(1 + A_1*A + A_2*G))*T*V - delta*T_i - psi_01*T_i*Z;
dy(3) = p*T_i - c*V - mu_v*G*V;
dy(4) = lambda_z + b*T_i*Z - mu_z*Z;
dy(5) = lambda_A - gamma_A*A*V - mu_A*A;
dy(6) = lambda_G + gamma_G*G*V - mu_G*G;

 
end