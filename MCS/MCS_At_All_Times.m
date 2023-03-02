clear all
close all
clc


numiter = 1000;


true_params = [47513.2875763767	0.170469182824029	5.94738085530006e-09	0.0315778580964400	6112.02230276937	8.10898357189916...
    0.200186087005715	4.04577582313757e-07	0.161802104330706	28.5877899533863	2.02925686877366e-08	8.74126737296193...
    0.0436946660173884	1.55765066550777e-12	0.0163710975483604	2.48330809208600e-09];

tforward = 0:0.1:170;

initial_cond = [1304000 0 200 10 3.44 2.57179227742938];


  
  
  [t, y_trp] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,true_params),...
      tforward,initial_cond);
  
  X = zeros(length(true_params),numiter);
  
 noiselevel = [0 0.01, 0.05, 0.1, 0.2];
 total_ARE =  zeros(length(noiselevel), length(true_params));
  
 total_ARE_Table = {'r', 'd', 'rho',  'delta', 'p', 'c', 'psi_01', 'b', ...
    'mu_z', 'lambda_A','gamma_A', 'mu_A', 'lambda_G', 'gamma_G', 'mu_G',...
    'mu_v'};
  
 save('Workspace_MCS_At_All_Times')

for j = 1: length(noiselevel)
    
rng default
noise = noiselevel(j);

    parfor i = 1:numiter
            
            fprintf('Current noise level is %g\n',noise);
            fprintf('Current iteration is %d\n',i);

    Viral_data = log10(y_trp(:,3));
    CD4_data =  log10(y_trp(:,1));
    Albumin_data = y_trp(:,5);
    Globulin_data = y_trp(:,6);

  ViralData = (noise*randn(length(tforward),1))+ Viral_data;
  CD4Data = (noise*randn(length(tforward),1))+ CD4_data;
  AlbuminData   = (noise*randn(length(tforward),1))+ Albumin_data; 
  GlobulinData   = (noise*randn(length(tforward),1))+ Globulin_data;
  
            k = true_params;
            lb = zeros(1,length(true_params));
             
             k = fminsearchbnd(@(k)err_in_data(k,ViralData, CD4Data,AlbuminData,GlobulinData),...
                 k,lb,[],optimset('MaxFunEvals', 1e+4,'MaxIter',1e+4));
       
             
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for p = 1:length(true_params)
        
        arescore(p) = 100*sum(abs(true_params(p) - X(p,:))/abs(true_params(p)))/numiter;
        
    end
    
    total_ARE(j,:) = round(arescore,1);
    total_ARE_Table(j+1,:) = num2cell(total_ARE(j,:));


end
  
  

function error_in_data = err_in_data(k,ViralData,CD4Data,AlbuminData,GlobulinData) 
 
tforward = 0:0.1:170;

initial_cond = [1304000 0 200 10 3.44 2.57179227742938];
 
  
  [~,y] = ode23s(@(t,y)Model_HIV_Control_WithinHost(y,k),tforward,initial_cond);
  

  
 Model_Viral = log10(y(:,3));
 Model_CD4 = log10(y(:,1));
 Model_Albumin = y(:,5); 
 Model_Globulin = y(:,6); 
 
 error_in_data = sum((Model_Viral - ViralData).^2) +...
              20*sum((Model_CD4 - CD4Data).^2)+...
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
