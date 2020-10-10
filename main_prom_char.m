close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis of extended  models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_names= {'N','E','A', 'I','R'}; 
u_names= {'m'}; % the input 
p_ext_names= {'k','Nm','d','kE','dE','vA','dA','kRIrel','k_basal','k_regulated','theta','eta'};

%XX define initial condition X0, time span tspan, and a vector m_vector for the changeable input 
x0 = [100000 0 0 0 0];
tspan = 62;
m_vector = (0:0.3:3);

%define reference parameter values from paper and introduce vectors for
%additional unknown parameters
p_ref([3 4 5 6])= [0.004 5 2 4.8*10^(-7)]; 
p_ref([1 2 7])= [0.970 1.24*10^9 0.639];
p_ref(8)= p_ref(5);
k_basal_vector = (0.12:0.05:0.32);
k_regulated_vector= (5:0.75:8);
theta_vector= (0.4:0.3:1.6);
eta_vector= (1.8:0.125:2.3);
choice = 3;

%two possible tasks here: selecting a characterization circuit or fitting
%the circuit's behavior to data and infering promoter parameter values 
%adapts to your needs
select_char_circuit= 0;
infer_param_values= 1;

if select_char_circuit
     for m_index = 1:length(m_vector)
%     for m_index = 1:2
        % for each inducer concentration: compute and store the variance of the output of the circuit for all possible parameter concentrations
        m = m_vector(m_index);
        output_youI = zeros(5,5,5,5); %initializing a 4D array to store all possible E values for different combinations, Storing 0 in all rn
        output_youR = zeros(5,5,5,5);
        output_youRI = zeros(5,5,5,5);
        run = 0; % just to keep a check on which run is going on since 625 permutations are possible for 4 set of parameters
        %XX here you iterate over the 4 parameter vectors, define p_ext as
        %appropriate, do the numerical simulation for the three circuits and 
        %store the appropriate value in 4-dimensional arrays output_youR,
        %output_youI and output_youRI.
        for kb = 1:length(k_basal_vector) %to do num sim for diff K_basal
            for kr = 1:length(k_regulated_vector) %to do num sim for diff K_regulated
                for theta = 1:length(theta_vector) %to do num sim for diff theta
                    for eta = 1:length(eta_vector) %to do num sim for diff eta
                        p = [p_ref(1:8) k_basal_vector(kb) k_regulated_vector(kr) theta_vector(theta) eta_vector(eta)];
                        Tend = tspan;
                        %if choice == 1
                        [t x_onI] = ode15s(@you_odeI, (0:Tend), x0, [], p,m);
                        output_youI(kb,kr,theta,eta) = x_onI(end,2); % to store values generated
                        %elseif choice==2
                        [t x_onR] = ode15s(@you_odeR, (0:Tend), x0, [], p,m);  
                        output_youR(kb,kr,theta,eta) = x_onR(end,2);
                        %elseif choice == 3
                        [t x_onRI] = ode15s(@you_odeRI, (0:Tend), x0, [], p,m);
                        output_youRI(kb,kr,theta,eta) = x_onRI(end,2);
                        %end
                        run = run+1; % again to display which run is going on
                        display(run);                        
                    end
                end
            end
        end
        %Note: if a is a multidimensional array (eg a= [1 2; 3 4]), 
        %then a(:) is a vector containing all the values in a 
        %(eg a(:)= [1 2 3 4])
        m_var_youI(m_index)= var(output_youI(:)); %variance of the outputs at a given inducer concentration and for different parameter values
        m_var_youR(m_index)= var(output_youR(:));
        m_var_youRI(m_index)= var(output_youRI(:));
    end
    var_youI= mean(m_var_youI); %overall variance defined as the mean over the variances for given inducer concentrations
    var_youR= mean(m_var_youR);
    var_youRI= mean(m_var_youRI);
    save('var_you.mat', 'var_youR', 'var_youI', 'var_youRI');
    display('Done');% to display that experiment is over
end

%here, take a decision and ask for experimental data

%here we assume you have gotten data in "data.mat"
if infer_param_values
    load('dataRI.mat','data');
    figure();
    plot(m_vector,data,'x'); hold on 
    p_search = [0.2 5 1 2]';
    p_fixed = p_ref(1:8)';
    Tend = tspan;
    opts=cmaes();
    
    opts.MaxFunEvals= 1500;  
    opts.LBounds = zeros(size(p_search)); % enforces positivity of the tested parameters
    [p_min, cost_min]= cmaes('compute_cost', p_search, p_search*0.6, opts, x0, Tend, m_vector, p_fixed, data', 3);
    p_opt= [p_fixed(1:8)' p_min'];
    E_ss = generate_ss_data(m_vector, tspan, x0, p_opt, 3);
    plot(m_vector,E_ss); hold on;
    xlabel('Concentration (m)')
    ylabel('E steady state')
    legend('Experimental', 'Simulated')
    title('Parameter fit for E_{ss} at different promoter concentrations')
    %save the optimized parameter values if you are happy with the search result. Comment line below otherwise
    save('p_opt.mat','p_opt');
end
   