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
p_ref([1 2 7])= [0.970 1.24*10^9 0.639];% set values for k, Nm and dA at pH 7
p_ref(8) = p_ref(5); % in the case of the ON system, kIRrel is equal to dE

p = p_ref; %simplifying so we can just call p instead of p_ref

%define initial condifions and simulation time (somewhat arbitrary)
tspan= 62; % like in the paper (in hrs);
x0= [100000 0 0 0 0]; % the initial condition of the system, assuming the preculture is made without IPTG (OFF case)

analyze_original_system= 1; % adapt based on your needs
analyze_modified_systems= 0;
use_default_params= 1; 
use_optimized_params= 0;
choice = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis of the original model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if analyze_original_system
    Tend = tspan;   
    [t x_on] = ode15s(@you_ode, (0:Tend), x0, [], p);
    %x_on(:,1);
    p_off = p;
    p_off(8) = 0; %turn system OFF
    [t x_off] = ode15s(@you_ode, (0:Tend), x0, [], p_off);
    %x_off(:,1)
    
    figure(1);
    semilogy(t,x_on(:,1),'-'); hold on
    semilogy(t,x_off(:,1),'-x'); hold off
    legend('N-on','N-off'); 
    xlabel('Time (Hrs)');
    ylabel('N (ml^-1)');
    title('Change in cell density with Time');
    
    figure(2);
    plot(t,x_on(:,2),'-',t,x_on(:,3),'.-'); hold on
    plot(t,x_off(:,2),'o',t,x_off(:,3),'-*'); hold off
    legend('E-on','A-on','E-off','A-off');
    xlabel('Time (Hrs)');
    ylabel('E and R (nM)');
    %yyaxis right
    %ylabel('A (nM)');
    title('Change in concentrations of A and E with Time');    
    
    figure(3);
    plot(t,x_on(:,4),'-',t,x_on(:,5),'.'); hold on
    plot(t,x_off(:,4),'o',t,x_off(:,5),'-*'); hold off
    legend('I-on','R-on','I-off','R-off');
    xlabel('Time (Hrs)');
    ylabel('Relative I and R');
    %yyaxis right
    %ylabel('R (nM)');
    title('Change in concentrations of I and R with Time');    
    %currentFigure = gcf;
    %title(currentFigure.Children(end), 'Analysis of the original circuit');
    
    p_PH62 = [0.885 1.25*10^(9) p(3:6) 0.274 p(8)];
    p_PH78 = [0.936 1.2*10^(9) p(3:6) 1.19 p(8)];
    Nmin = ((p_PH62(1)*p_PH62(5)*p_PH62(7))/(p_PH62(6)*p_PH62(4)*p_PH62(3))); % new Ns formula in accordance with new ODEs
    Nmax = ((p_PH78(1)*p_PH78(5)*p_PH78(7))/(p_PH78(6)*p_PH78(4)*p_PH78(3)));
    dynamic_range = Nmax/Nmin;
    display(dynamic_range);
    
    % to plot system dynamics
    p_PH66 = [0.928 1.17*10^(9) p(3:6) 0.304 p(8)];
    p_PH74 = [0.897 1.16*10^(9) p(3:6) 0.791 p(8)];
    [t x_on1] = ode15s(@you_ode, (0:Tend), x0, [], p_PH62);
    [t x_on2] = ode15s(@you_ode, (0:Tend), x0, [], p_PH66);
    [t x_on3] = ode15s(@you_ode, (0:Tend), x0, [], p);
    [t x_on4] = ode15s(@you_ode, (0:Tend), x0, [], p_PH74);
    [t x_on5] = ode15s(@you_ode, (0:Tend), x0, [], p_PH78);
    figure(4);
    plot(t,x_off(:,1));hold on
    plot(t,x_on1(:,1));
    plot(t,x_on2(:,1));
    plot(t,x_on3(:,1));
    plot(t,x_on4(:,1));
    plot(t,x_on5(:,1)); 
    legend('OFF','6.6','6.8','7','7.4','7.8');
    ylabel('N (CFU-/ml) at different pH')
    xlabel('Time (hrs)')
    title('Behaviour of the system at different pH')
end

if analyze_modified_systems
    %extended parameter set
    p_ext_names= {'k','Nm','d','ke','de','va','da','kRIrel','k_basal','k_regulated','theta','eta'};
    if use_default_params
        p_ref([3 4 5 6])= [0.004 5 2 4.8*10^(-7)]; %set values for d, kE, dE, and vA
        p_ref([1 2 7])= [0.970 1.24*10^9 0.639];% set values for k, Nm and dA 
        p_ref(8) = p_ref(5);
        p_default(1:8)= p_ref;
        p_default(9)= 0.2; %k_basal
        p_default(10)= 5; %k_regulated
        p_default(11)= 1; %theta
        p_default(12)= 2; %eta
        p= p_default;
    elseif use_optimized_params
        load('p_opt.mat', 'p_opt');
        p = p_opt;
    end
    
    m_vector = (0:0.3:3); %the two extreme values of m with increments of 0.3
    if choice == 1
    %do numerical simulation for the "I" circuit
    %Here do num sim for the ON case with m absent and store in x_on 
    Tend = tspan;
    m = m_vector(1); % m = 0 if m_vector(1)
    [t x_on] = ode15s(@you_odeI, (0:Tend), x0, [], p, m);
    Nmax = x_on(end,1);
    
    %XX here do num sim for the OFF case with m absent and store in x_off 
    p_off = p;
    p_off(8) = 0;
    m = m_vector(1);
    [t x_off] = ode15s(@you_odeI, (0:Tend), x0, [], p_off,m);
    x_off(end,1);
    
    %figures
    figure(1);
    semilogy(t,x_on(:,1),'-'); hold on
    semilogy(t,x_off(:,1),'-x'); hold off
    legend('N-on','N-off'); 
    xlabel('Time (Hrs)');
    ylabel('N (ml^-1)');
    title('Analysis of I Circuit: Change in cell density with Time at m=0');
    
    figure(2);
    plot(t,x_on(:,2),'-',t,x_on(:,3),'-x'); hold on
    plot(t,x_off(:,2),'o',t,x_off(:,3),'*'); hold off
    legend('E-on','A-on','E-off','A-off');
    xlabel('Time (Hrs)');
    ylabel('E (nM)');
    yyaxis right
    ylabel('A (nM)');
    title('Change in concentrations of A and E with Time at m=0');   
    
    figure(3);
    plot(t,x_on(:,4),'-',t,x_on(:,5),'-x'); hold on
    plot(t,x_off(:,4),'o',t,x_off(:,5),'*'); hold off
    legend('I-on','R-on','I-off','R-off');    
    xlabel('Time (Hrs)');
    ylabel('I relative');
    yyaxis right
    ylabel('R relative');
    title('Change in concentrations of I and R with Time at m=0');    
    %title(currentFigure.Children(end), 'Analysis of the modified I circuit');
    
    %figure to compare different curves at diff (m)
    figure(4);
    for varying_m=(0:0.3:3)
        [t x_on] = ode15s(@you_odeI, (0:Tend), x0, [], p, varying_m);
        semilogy(t, x_on(:,1), '-'); hold on;
        %display(varying_m)
        %display(x_on(end, 1))
        legend('m=0.0','m=0.3','m=0.6','m=0.9','m=1.2','m=1.5','m=1.8','m=2.1','m=2.4','m=2.7','m=3.0');
        xlabel('Time (Hrs)');
        ylabel('N (ml^-1)');
        title('Analysis of I Circuit: Change in Cell Density with Time and Varying m');
    end
        
    if use_default_params
        p= p_default;
    elseif use_optimized_params
        p= p_opt;
    end
    
    %XX here do num sim for the ON case with m present in large concentration and store in x_on 
    Tend = tspan;
    m = m_vector(end); % m = 0 if m_vector(1)
    [t x_on1] = ode15s(@you_odeI, (0:Tend), x0, [], p, m);
    
    %XX here do num sim for the OFF case with m present in large concentration and store in x_off 
    p_off = p;
    p_off(8) = 0;
    m = m_vector(end);
    [t x_off1] = ode15s(@you_odeI, (0:Tend), x0, [], p_off,m);
    x_off1(end,1);
    
    %XX here define Nmin or Nmax, depending on what is appropriate
    Nmin = x_on1(end,1);

    %figures
    figure(5);
    semilogy(t,x_on1(:,1),'-'); hold on
    semilogy(t,x_off1(:,1),'-x'); hold off
    legend('N-on','N-off'); 
    xlabel('Time (Hrs)');
    ylabel('N (ml^-1)');
    title('Analysis of I Circuit: Change in cell density with Time at m=3');
    
    figure(6);
    plot(t,x_on1(:,2),'-',t,x_on(:,3),'-x'); hold on
    plot(t,x_off1(:,2),'o',t,x_off1(:,3),'*'); hold off
    legend('E-on','A-on','E-off','A-off');
    xlabel('Time (Hrs)');
    ylabel('E (nM)');
    yyaxis right
    ylabel('A (nM)');
    title('Change in concentrations of A and E with Time at m=3');   
    
    figure(7);
    plot(t,x_on1(:,4),'-',t,x_on1(:,5),'-x'); hold on
    plot(t,x_off1(:,4),'o',t,x_off1(:,5),'*'); hold off
    legend('I-on','R-on','I-off','R-off'); 
    xlabel('Time (Hrs)');
    ylabel('I relative');
    yyaxis right
    ylabel('R relative');
    title('Change in concentrations of I and R with Time at m=3');    
    
    %define and display the fold change 
    dynamic_range= Nmax/Nmin;
    display(dynamic_range);
    
    % XX do similar tests with the "R" and "RI" circuits
    elseif choice == 2
    %do numerical simulation for the "R" circuit
    %XX here do num sim for the ON case with m absent and store in x_on 
    Tend = tspan;
    m = m_vector(1); % m = 0 if m_vector(1)
    [t x_on2] = ode15s(@you_odeR, (0:Tend), x0, [], p, m);
    Nmax1 = x_on2(end,1);
    
    %XX here do num sim for the OFF case  with m absent and store in x_off 
    p_off = p;
    p_off(8) = 0;
    m = m_vector(1);
    [t x_off2] = ode15s(@you_odeR, (0:Tend), x0, [], p_off,m);
    x_off2(end,1);
    
    %figures
    figure(1);
    semilogy(t,x_on2(:,1),'-'); hold on
    semilogy(t,x_off2(:,1),'-x'); hold off
    legend('N-on','N-off'); 
    xlabel('Time (Hrs)');
    ylabel(['N (ml^-1)']);
    title('Analysis of R Circuit: Change in cell density with Time at m=0');
    figure(2);
    plot(t,x_on2(:,2),'-',t,x_on2(:,3),'-x'); hold on
    plot(t,x_off2(:,2),'o',t,x_off2(:,3),'*'); hold off
    legend('E-on','A-on','E-off','A-off');
    title('Change in concentrations of A and E with Time at m=0');   

    figure(3);
    plot(t,x_on2(:,4),'-',t,x_on2(:,5),'.-'); hold on
    plot(t,x_off2(:,4),'o',t,x_off2(:,5),'*'); hold off
    legend('I-on','R-on','I-off','R-off');    
    title('Change in concentrations of I and R with Time at m=0');   
       
    if use_default_params
        p= p_default;
    elseif use_optimized_params
        p= p_opt;
    end
    %XX here do num sim for the ON case with m present in large concentration and store in x_on 
    Tend = tspan;
    m = m_vector(end); % m = 0 if m_vector(1)
    [t x_on3] = ode15s(@you_odeR, (0:Tend), x0, [], p, m);
    
    %XX here do num sim for the OFF case with m present in large concentration and store in x_off 
    p_off = p;
    p_off(8) = 0;
    m = m_vector(end);
    [t x_off3] = ode15s(@you_odeR, (0:Tend), x0, [], p_off,m);
    x_off3(end,1);
    
    %XX here define Nmin or Nmax, depending on what is appropriate
    Nmin1 = x_on3(end,1);

    %XX here is a figure block [update the title!]
    %figures
    figure(4);
    semilogy(t,x_on3(:,1),'-'); hold on
    semilogy(t,x_off3(:,1),'-x'); hold off
    legend('N-on','N-off'); 
    xlabel('Time (Hrs)');
    ylabel('N (ml^-1)');
    title('Analysis of R Circuit: Change in cell density with Time at m=3');
    
    figure(5);
    plot(t,x_on3(:,2),'-',t,x_on3(:,3),'-x'); hold on
    plot(t,x_off3(:,2),'o',t,x_off3(:,3),'*'); hold off
    legend('E-on','A-on','E-off','A-off');
    xlabel('Time (Hrs)');
    ylabel('E (nM)');
    yyaxis right
    ylabel('A (nM)');
    title('Change in concentrations of A and E with Time at m=3'); 
    
    figure(6);
    plot(t,x_on3(:,4),'-',t,x_on3(:,5),'x-'); hold on
    plot(t,x_off3(:,4),'o',t,x_off3(:,5),'*'); hold off
    legend('I-on','R-on','I-off','R-off'); 
    xlabel('Time (Hrs)');
    ylabel('I (nM)');
    yyaxis right
    ylabel('R (nM)');
    title('Change in concentrations of I and R with Time at m=3');    
    
    %define and display the fold change 
    dynamic_range1= Nmax1/Nmin1;
    display(dynamic_range1);
    
    elseif choice == 3
    %XX IR circuit
    %XX here do num sim for the ON case with m absent and store in x_on 
    Tend = tspan;
    m = m_vector(1); % m = 0 if m_vector(1)
    [t x_on4] = ode15s(@you_odeRI, (0:Tend), x0, [], p, m);
    Nmax2 = x_on4(end,1); 
    
    %XX here do num sim for the OFF case  with m absent and store in x_off 
    p_off = p;
    p_off(8) = 0;
    m = m_vector(1);
    [t x_off4] = ode15s(@you_odeRI, (0:Tend), x0, [], p_off,m);
    x_off4(end,1);
    
    %figures
    figure(1);
    semilogy(t,x_on4(:,1),'-'); hold on
    semilogy(t,x_off4(:,1),'-x'); hold off
    legend('N-on','N-off');  
    xlabel('Time (Hrs)');
    ylabel('N (ml^-1)');
    title('Analysis of RI Circuit: Change in cell density with Time at m=0');
    
    figure(2);
    plot(t,x_on4(:,2),'-',t,x_on4(:,3),'x-'); hold on
    plot(t,x_off4(:,2),'o',t,x_off4(:,3),'*'); hold off
    legend('E-on','A-on','E-off','A-off');
    xlabel('Time (Hrs)');
    ylabel('E (nM)');
    yyaxis right
    ylabel('A (nM)');
    title('Change in concentrations of A and E with Time at m=0');   
    
    figure(3);
    plot(t,x_on4(:,4),'-',t,x_on4(:,5),'.-'); hold on
    plot(t,x_off4(:,4),'o',t,x_off4(:,5),'*'); hold off
    legend('I-on','R-on','I-off','R-off');    
    xlabel('Time (Hrs)');
    ylabel('I (nM)');
    yyaxis right
    ylabel('R (nM)');
    title('Change in concentrations of I and R with Time at m=0');    
   
       
    if use_default_params
        p= p_default;
    elseif use_optimized_params
        p= p_opt;
    end
    %XX here do num sim for the ON case with m present in large concentration and store in x_on 
    Tend = tspan;
    m = m_vector(end); % m = 0 if m_vector(1)
    [t x_on5] = ode15s(@you_odeRI, (0:Tend), x0, [], p, m);
    
    %XX here do num sim for the OFF case with m present in large concentration and store in x_off 
    p_off = p;
    p_off(8) = 0;
    m = m_vector(end);
    [t x_off5] = ode15s(@you_odeRI, (0:Tend), x0, [], p_off,m);
    x_off5(end,1);
    
    %XX here define Nmin or Nmax, depending on what is appropriate
    Nmin2 = x_on5(end,1);

    %figures
    figure(4);
    semilogy(t,x_on5(:,1),'-'); hold on
    semilogy(t,x_off5(:,1),'-x'); hold off
    legend('N-on','N-off');  
    xlabel('Time (Hrs)');
    ylabel('N (ml^-1)');
    title('Analysis of RI Circuit: Change in cell density with Time at m=3');
    
    figure(5);
    plot(t,x_on5(:,2),'-',t,x_on5(:,3),'-x'); hold on
    plot(t,x_off5(:,2),'o',t,x_off5(:,3),'*'); hold off
    legend('E-on','A-on','E-off','A-off');
    xlabel('Time (Hrs)');
    ylabel('E (nM)');
    yyaxis right
    ylabel('A (nM)');
    title('Change in concentrations of A and E with Time at m=3'); 
    
    figure(6);
    plot(t,x_on5(:,4),'-',t,x_on5(:,5),'.-'); hold on
    plot(t,x_off5(:,4),'o',t,x_off5(:,5),'*'); hold off
    legend('I-on','R-on','I-off','R-off'); 
    xlabel('Time (Hrs)');
    ylabel('I (nM)');
    yyaxis right
    ylabel('R (nM)');
    title('Change in concentrations of I and R with Time at m=3');    
    
    %define and display the fold change 
    dynamic_range2= Nmax2/Nmin2;
    display(dynamic_range2);  
    end
end    

       