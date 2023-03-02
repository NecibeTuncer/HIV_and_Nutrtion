clear all
close all
clc

global  tforward initial_cond ViralData  CD4Data AlbuminData GlobulinData

ViralData = [8.51E+6	1.13E+7	9.17E+5	3.82E+5	2.73E+5	1.57E+5	1.55E+5	9.51E+4]';


CD4Data = 1000*[1304	619.5	371.1666667	151.6666667	235.8333333	241.1666667	270	379.1666667	246.1666667]';

AlbuminData = [3.44299233423559;2.94258941669938;3.12718763559725;...
               2.92043928342666;3.87677955245160;3.22351592040829;...
               3.51332727235159;3.34054796784925;3.06362997706469]; %mg/dL

GlobulinData = [2.57179227742938;2.59199354586829;2.70640482476377];

initial_cond = [1304000 0 200 10 3.44 2.57179227742938];
tVLdata = [7,2*7,3*7,4*7,6*7,8*7,10*7,12*7];
tCD4data = [0,7,2*7,3*7,4*7,6*7,8*7,10*7,12*7];
tAlbumindata = [0, 7, 7*2, 7*3,5*7,6*7, 8*7, 12*7, 16*7];  
tGlobulindata = [0, 12*7, 24*7];
tforward = 0:0.1:170;




k = [38732.0471628024	0.146678171205492	5.41871610765662e-09	1.02575918493109	5919.22197096820	8.08364564189790...
    0.0617916954117044	3.87301386642762e-07	0.0646390842699867	48.3029843480593	6.65025922423098e-08	14.6899133017896...
    0.0235274962694968	1.72786049371802e-12	0.00870379317329876	1.02068289821756e-09];


%results in 

% k =[47513.2875763767	0.170469182824029	5.94738085530006e-09	0.0315778580964400	6112.02230276937	8.10898357189916...
%     0.200186087005715	4.04577582313757e-07	0.161802104330706	28.5877899533863	2.02925686877366e-08	8.74126737296193...
%     0.0436946660173884	1.55765066550777e-12	0.0163710975483604	2.48330809208600e-09];

lb = zeros(1,length(k));



[k,fval] =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter','MaxFunEvals', 1e+4, 'MaxIter', 1e+4));


[t_r, y_r] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,k),tforward,initial_cond);

display('Parameters after data fitting Model 1:')

fprintf('r = %g\n',  k(1));   
fprintf('d = %g\n', k(2));
fprintf('rho = %g\n', k(3));
fprintf('delta = %g\n',  k(4));
fprintf('p = %g\n', k(5));
fprintf('c = %g\n',  k(6));
fprintf('psi_01 = %g\n',  k(7));
fprintf('b = %g\n',  k(8));
fprintf('mu_z = %g\n',  k(9));
fprintf('lambda_A = %g\n', k(10));
fprintf('gamma_A = %g\n', k(11));
fprintf('mu_A = %g\n', k(12));
fprintf('lambda_G = %g\n', k(13));
fprintf('gamma_G = %g\n', k(14));
fprintf('mu_G = %g\n', k(15));
fprintf('mu_v = %g\n', k(16));
%fprintf('lambda_z = %g\n', k(17));
%fprintf('SSE = %g\n', fval);




figure(1)
plot(tforward, log10(y_r(:,3)),'-b','LineWidth',4,'color', 'b')
hold on 
plot(tVLdata, log10(ViralData), '.', 'color', 'r', 'MarkerSize', 40)
title('HIV Viral Load' , 'vRNA copies/mL', 'fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('Log_{10} V(t)','fontweight', 'normal', 'fontsize', 18)
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3, 'fontsize',18)
xlim([0,175])
yticks([2 4 6 8])
yticklabels({'10^2','10^4','10^6', '10^8'})

figure(5)
semilogy(tforward, y_r(:,3),'-b','LineWidth',4,'color', 'b')
hold on 
plot(tVLdata, ViralData, '.', 'color', 'r', 'MarkerSize', 40)
title('HIV Viral Load' , 'vRNA copies/mL', 'fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('Log_{10} V(t)','fontweight', 'normal', 'fontsize', 18)
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 1, 'fontsize',18)
xlim([0,175])

figure(2)
plot(tforward, log10(y_r(:,1)),'-b','LineWidth',4, 'color', 'b')
hold on 
plot(tCD4data, log10(CD4Data), '.', 'color', 'r', 'MarkerSize', 40)
title('CD4 Cell Count','CD4 cells/mL', 'fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('Log_{10} T(t)', 'fontweight', 'normal', 'fontsize', 18)
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3,'fontsize',18)
xlim([0,175])

figure(3)
plot(tforward, y_r(:,5),'-b','LineWidth',4, 'color', 'b')
hold on
plot(tAlbumindata, AlbuminData, '.', 'color', 'r', 'MarkerSize', 40)
title('Albumin', 'mg/dL','fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('A(t)', 'fontweight', 'normal', 'fontsize', 18)
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3,'fontsize',18)
xlim([0,175])

figure(4)
plot(tforward, y_r(:,6),'-b','LineWidth',4, 'color', 'b')
hold on
plot(tGlobulindata, GlobulinData, '.', 'color', 'r', 'MarkerSize', 40)
title('Globulin', 'mg/dL' , 'fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('G(t)', 'fontweight', 'normal', 'fontsize', 18) %what are units?
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',18)
xlim([0,175])

 function error_in_data = err_in_data(k) 
 
 global  tforward initial_cond ViralData CD4Data  AlbuminData GlobulinData
 
  
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
rho = k(3);
delta = k(4);
p = k(5);
c = k(6);

psi_01 = k(7);
b = k(8);
mu_z = k(9);

lambda_A = k(10); 
gamma_A = k(11);
mu_A = k(12); 

lambda_G = k(13); 
gamma_G =  k(14);
mu_G = k(15); 
mu_v = k(16);

lambda_z = 1;




T = y(1);
T_i = y(2);
V = y(3);
Z = y(4);
A = y(5);
G = y(6);




dy(1) = r - rho*T*V - d*T;
dy(2) = rho*T*V - delta*T_i - psi_01*T_i*Z;
dy(3) = p*T_i - c*V - mu_v*G*V;
dy(4) = lambda_z + b*T_i*Z - mu_z*Z;
dy(5) = lambda_A - gamma_A*A*V - mu_A*A;
dy(6) = lambda_G + gamma_G*G*V - mu_G*G;

 
end