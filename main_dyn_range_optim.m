close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis of the original model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_names= {'N','E','A','I','R'}; % N is in CFU
p_names= {'k','Nm','d','kE','dE','vA','dA','kRIrel'};
%p(8), called here kRIrel is the protein production rate parameter for I
%and R proteins. It is such that at steady state in the ON case, the level of the I and R
%proteins is 1, given that the degradation rate parameter of these proteins
%is the same as for the E protein, dE;

%define reference parameter values (from paper)
p_ref([3 4 5 6])= [0.004 5 2 4.8*10^(-7)]; %set values for d, kE, dE, and vA
p_ref([1 2 7])= [0.970 1.24*10^9 0.639]; %set values for k, Nm and dA at pH 7
p_ref(8) = p_ref(5); %in the case of the ON system, kIRrel is equal to dE

p = p_ref; %simplifying so we can just call p instead of p_ref

%define initial conditions and simulation time (somewhat arbitrary)
tspan= 62; % they use 62 hours in the paper, use 80 to extend system time 
x0= [100000 0 0 0 0]; % the initial condition of the system, assuming the preculture is made without IPTG (OFF case)

analyze_original_system= 1; % adapt based on your needs
analyze_modified_systems= 0;
use_default_params= 0; 
use_optimized_params= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis of the original model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analyze_original_system
    % Computes ON and OFF behaviors of original system; plot on a same figure
    % % 
    %XX set parameters p for the ON system (adapts parameters to mimic the presence of IPTG) and do numerical simulation
    Tend= tspan;   
    [t x_on] = ode15s(@you_ode, [0, Tend], x0, [], p);
    %XX set parameters p for the OFF system and do numerical simulation   
    p_off = p;
    p_off(8) = 0; 
    [t x_off] = ode15s(@you_ode, [0, Tend], x0, [], p_off);
    figure();
    %I call the 8 lines below a figure block (reused later)
    subplot(3,1,1); semilogy(t,x_on(:,1),'-'); legend('N-on'); hold on 
    subplot(3,1,2); plot(t,x_on(:,2:3),'-'); legend('E-on','A-on'); hold on
    subplot(3,1,3); plot(t,x_on(:,4:5),'-'); legend('I-on','R-on'); hold on
    subplot(3,1,1); semilogy(t,x_off(:,1),'--'); legend('N-off'); 
    subplot(3,1,2); plot(t,x_off(:,2:3),'--'); legend('E-off','A-off'); 
    subplot(3,1,3); plot(t,x_off(:,4:5),'--'); legend('I-off','R-off'); hold on
    currentFigure = gcf;
    %title(currentFigure.Children(end), 'Analysis of the original circuit');
    
    % Compute cell density of engineered system at pH6.2 and pH7.8 and the corresponding dynamic range
    
    %p_pH62=XX;% stores values of k, Nm and da at pH 6.2   
    p_PH62 = [0.885 1.25*10^(9) p(3:6) 0.274 0];
    %p_pH78=XX;% same at pH 7.8   
    p_PH78 = [0.936 1.2*10^(9) p(3:6) 1.19 0];
    %XX compute cell density at steady state for pH6.2; stored in Nmin
    Nmin = (p_PH62(1)*p_PH62(5)*p_PH62(7))/(p_PH62(6)*p_PH62(4)*p_PH62(3));
    %XX compute cell density at steady state for pH7.8; stored in Nmax
    Nmax = (p_PH78(1)*p_PH78(5)*p_PH78(7))/(p_PH78(6)*p_PH78(4)*p_PH78(3));
    %define and display the fold change associated to pH variations
    dynamic_range = Nmax / Nmin;
    display(dynamic_range);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis of extended  models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analyze_modified_systems
    %extended parameter set
    p_ext_names= {'k','Nm','d','ke','de','va','da','kRIrel','k_basal','k_regulated','theta','eta'};
    if use_default_params
        p_ref([3 4 5 6])= XX; %set values for d, kE, dE, and vA to pH 7?
        p_ref([1 2 7])= XX;% set values for k, Nm and dA 
        p_ref(8)= XX;
        p_default(1:8)= XX;
        p_default(9)= 0.2; %k_basal
        p_default(10)= 5; %k_regulated
        p_default(11)= 1; %theta
        p_default(12)= 2; %eta
        p= p_default;
    elseif use_optimized_params
        load('p_opt.mat', 'p_opt');
        p= p_opt;
    end
    m_vector= XX; %the two extreme values of m
    %do numerical simulation for the "I" circuit
    %XX here do num sim for the ON case with m absent and store in x_on 
    %XX here do num sim for the OFF case  with m absent and store in x_off 
    %XX here define Nmin or Nmax, depending on what is appropriate
    figure();
    %XX here is a figure block [update the title!]
    if use_default_params
        p= XX
    elseif use_optimized_params
        p= XX
    end
    %XX here do num sim for the ON case with m present in large concentration and store in x_on 
    %XX here do num sim for the OFF case with m present in large concentration and store in x_off 
    %XX here define Nmin or Nmax, depending on what is appropriate
    figure();
    %XX here is a figure block [update the title!]
    %define and display the fold change 
    dynamic_range= XX;
    display(dynamic_range);
    % XX do similar tests with the "R" and "RI" circuits
    % that is copy/paste and adapt the above block for the other two
    % circuits
end    
   