%prediction_plot
clear all
close all
clc


% params = [47513.2875763767	0.170469182824029	5.94738085530006e-09	0.0315778580964400	6112.02230276937	8.10898357189916...
%     0.200186087005715	4.04577582313757e-07	0.161802104330706	28.5877899533863	2.02925686877366e-08	8.74126737296193...
%     0.0436946660173884	1.55765066550777e-12	0.0163710975483604	2.48330809208600e-09];


lambda_A = 3 + rand(1,100)*(300-3);
lambda_G = 0.005 + rand(1,100)*(0.5 - 0.005);

n = length(lambda_A);

initial_cond = [1304000 0 200 10 3.44 2.57179227742938];
tforward = 0:0.1:170;

%quartile_limits =  [0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995];

quartile_limits =  [0.05,0.25,0.5,0.75,0.95];

m = length(quartile_limits);

Viral_predictions = zeros(n,length(tforward));
Albumin_predictions = zeros(n,length(tforward));
Globulin_predictions = zeros(n,length(tforward));

  for i = 1:n
      
      params = [lambda_A(i) lambda_G(i)] ;
      
      [~, y_trp] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,params),...
      tforward,initial_cond);
   
      
      Viral_predictions(i,:) =  log10(y_trp(:,3));
      Albumin_predictions(i,:) =  y_trp(:,5);
      Globulin_predictions(i,:) =  y_trp(:,6);

  end
  
  clear i
  clear j
  
    for i = 1:m
  
  p = quartile_limits(i);
  
        for j = 1:length(tforward)
          
            
            Viral_prediction_quartiles(i,j) = interp1(sort(Viral_predictions(:,j)),(n-1)*p+1);
            Albumin_prediction_quartiles(i,j) = interp1(sort(Albumin_predictions(:,j)),(n-1)*p+1);
            Globulin_prediction_quartiles(i,j) = interp1(sort(Globulin_predictions(:,j)),(n-1)*p+1);
             
        end
    end


params = [28.5877899533863  0.0436946660173884];
params = [303/2  0.505/2];
   [~, y_mean] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,params),...
      tforward,initial_cond);
   
 


  dimc = [0.9 0.9 0.9];


 figure(1); 
 ax1 = axes('Position',[0.2 0.1 0.7 0.8]);  % xlocation, ylocation, xsize, ysize
 ax2 = axes('Position',[0.65 0.65 0.2 0.2]);

  fillyyax(ax1,tforward,Viral_prediction_quartiles(1,:),Viral_prediction_quartiles(5,:),dimc);
  hold on 
  fillyyax(ax1,tforward,Viral_prediction_quartiles(2,:),Viral_prediction_quartiles(4,:),dimc.*0.9);
  hold on
  hold on
  plot(ax1,tforward,log10(y_mean(:,3)),'k','LineWidth',2)
  annotation('rectangle',[0.23+.5/3 1.5/3 .05 .05])
  annotation('arrow',[0.29+.5/3  0.35+.5/3],[1.6/3 1.8/3]) %xrange yrange
  set(ax1,'FontSize',12,'FontName','Arial','linewidth',2,'FontWeight','normal')
  set(ax1, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
  xlabel(ax1,'Time (days)','FontSize',18,'FontName','Arial','FontWeight','normal');
  ylabel(ax1,'log_{10}V','FontSize',18,'FontName','Arial','FontWeight','normal');

  fillyyax(ax2,tforward,Viral_prediction_quartiles(1,:),Viral_prediction_quartiles(5,:),dimc);
  hold on 
  fillyyax(ax2,tforward,Viral_prediction_quartiles(2,:),Viral_prediction_quartiles(4,:),dimc.*0.9);
  hold on
  hold on
  plot(ax2,tforward,log10(y_mean(:,3)),'k','LineWidth',2)
  set(ax2,'FontSize',8,'FontName','Arial','linewidth',2,'FontWeight','normal','box','on')
  xlim(ax2, [60.1,60.12])
  ylim(ax2, [5.18822,5.1885])
  hold off
  


  figure(2); 
  fillyy(tforward,Albumin_prediction_quartiles(1,:),Albumin_prediction_quartiles(5,:),dimc);
  hold on 
  fillyy(tforward,Albumin_prediction_quartiles(2,:),Albumin_prediction_quartiles(4,:),dimc.*0.9);
  hold on
  %plot(data.l_time,data.sputum,'ro')
  hold on
  plot(tforward,y_mean(:,5),'k','LineWidth',2)
  set(gca,'FontSize',12,'FontName','Arial','linewidth',2,'FontWeight','normal')
  set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
  xlabel('Time (days)','FontSize',18,'FontName','Arial','FontWeight','normal');
  ylabel('A(t)','FontSize',18,'FontName','Arial','FontWeight','normal');
  hold off




    figure(3); 
  fillyy(tforward,Globulin_prediction_quartiles(1,:),Globulin_prediction_quartiles(5,:),dimc);
  hold on 
  fillyy(tforward,Globulin_prediction_quartiles(2,:),Globulin_prediction_quartiles(4,:),dimc.*0.9);
  hold on
  %plot(data.l_time,data.sputum,'ro')
  hold on
  plot(tforward,y_mean(:,6),'k','LineWidth',2)
  set(gca,'FontSize',12,'FontName','Arial','linewidth',2,'FontWeight','normal')
  set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
  xlabel('Time (days)','FontSize',18,'FontName','Arial','FontWeight','normal');
  ylabel('G(t)','FontSize',18,'FontName','Arial','FontWeight','normal');
  hold off

function dy = Model_HIV_Control_WithinHost(y,params)

dy = zeros(6,1);

k = [47513.2875763767	0.170469182824029	5.94738085530006e-09	0.0315778580964400	6112.02230276937	8.10898357189916...
    0.200186087005715	4.04577582313757e-07	0.161802104330706		2.02925686877366e-08	8.74126737296193...
    1.55765066550777e-12	0.0163710975483604	2.48330809208600e-09];

r = k(1);
d = k(2);
rho = k(3);
delta = k(4);
p = k(5);
c = k(6);

psi_01 = k(7);
b = k(8);
mu_z = k(9);

lambda_A = params(1); %k(10); 
gamma_A = k(10);
mu_A = k(11); 

lambda_G = params(2); %k(12); 
gamma_G =  k(12);
mu_G = k(13); 
mu_v = k(14);

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



function out = fillyy(x,y1,y2,col)




    if nargin < 4
        col='red';
    end

x  = x(:)';
y1 = y1(:)';
y2 = y2(:)';
n   = length(x);
X = [ x(1),  x,  x(n),  fliplr(x)  ];
Y = [ y1(1), y2, y1(n), fliplr(y1) ];
h = fill(X,Y,col,'Linestyle','none');

    if nargout>0
        out=h;
    end
end


function out = fillyyax(ax,x,y1,y2,col)

x  = x(:)';
y1 = y1(:)';
y2 = y2(:)';
n   = length(x);
X = [ x(1),  x,  x(n),  fliplr(x)  ];
Y = [ y1(1), y2, y1(n), fliplr(y1) ];
h = fill(ax,X,Y,col,'Linestyle','none');

  
end