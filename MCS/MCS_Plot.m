clear all
close all
clc


numiter = 1;



true_params = [47513.2875763767	0.170469182824029	5.94738085530006e-09	0.0315778580964400	6112.02230276937	8.10898357189916...
    0.200186087005715	4.04577582313757e-07	0.161802104330706	28.5877899533863	2.02925686877366e-08	8.74126737296193...
    0.0436946660173884	1.55765066550777e-12	0.0163710975483604	2.48330809208600e-09];

tforward = 0:0.1:170;

initial_cond = [1304000 0 200 10 3.44 2.57179227742938];

 t_v_measure = [7,2*7,3*7,4*7,6*7,8*7,10*7,12*7]./0.1 +1;
 t_cd4_measure = [0,7,2*7,3*7,4*7,6*7,8*7,10*7,12*7]./0.1 +1;
 t_alb_measure = [0, 7, 7*2, 7*3,5*7,6*7, 8*7, 12*7, 16*7]./0.1 +1;
 t_glob_measure = [0,12*7, 24*7]./0.1 +1;
  
 tVLdata = [7,2*7,3*7,4*7,6*7,8*7,10*7,12*7];
tCD4data = [0,7,2*7,3*7,4*7,6*7,8*7,10*7,12*7];
tAlbumindata = [0, 7, 7*2, 7*3,5*7,6*7, 8*7, 12*7, 16*7];  
tGlobulindata = [0, 12*7, 24*7];
  
  [t, y_r] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,true_params),...
      tforward,initial_cond);
  

  

noise = 0.0;

  
         

    ViralData = log10(y_r(t_v_measure(:),3));
    CD4_data =  log10(y_r(t_cd4_measure(:),1));
    Albumin_data = y_r(t_alb_measure(:),5);
    Globulin_data = y_r(t_glob_measure(:),6);

  ViralData = (noise*randn(length(t_v_measure),1))+ ViralData;
  CD4Data = (noise*randn(length(t_cd4_measure),1))+ CD4_data;
  AlbuminData   = (noise*randn(length(t_alb_measure),1))+ Albumin_data; 
  GlobulinData   = (noise*randn(length(t_glob_measure),1))+ Globulin_data;
  
  figure(1)
plot(tforward, log10(y_r(:,3)),'-b','LineWidth',3,'color', 'b')
hold on 
plot(tVLdata, ViralData, '.', 'color', 'r', 'MarkerSize', 35)
title('Viral Load' , 'vRNA copies/mL', 'fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('Log_{10}V(t)','fontweight', 'normal', 'fontsize', 18)
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3, 'fontsize',18)
xlim([0,175])

figure(2)
plot(tforward, log10(y_r(:,1)),'-b','LineWidth',3, 'color', 'b')
hold on 
plot(tCD4data, CD4Data, '.', 'color', 'r', 'MarkerSize', 35)
title('CD4 Cell Count','CD4 cells/mL', 'fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('Log_{10}T(t)', 'fontweight', 'normal', 'fontsize', 18)
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3,'fontsize',18)
xlim([0,175])

figure(3)
plot(tforward, y_r(:,5),'-b','LineWidth',3, 'color', 'b')
hold on
plot(tAlbumindata, AlbuminData, '.', 'color', 'r', 'MarkerSize', 35)
title('Albumin', 'mg/dL','fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('A(t)', 'fontweight', 'normal', 'fontsize', 18)
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3,'fontsize',18)
xlim([0,175])

figure(4)
plot(tforward, y_r(:,6),'-b','LineWidth',3, 'color', 'b')
hold on
plot(tGlobulindata, GlobulinData, '.', 'color', 'r', 'MarkerSize', 35)
title('Globulin', 'mg/dL' , 'fontweight', 'normal', 'fontsize', 18)
xlabel('Time in Days', 'fontweight', 'normal', 'fontsize', 18)
ylabel('G(t)', 'fontweight', 'normal', 'fontsize', 18) %what are units?
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',18)
xlim([0,175])




  
  

   
 
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
