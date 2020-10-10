function cost= compute_cost(p_search, x0, tspan, m_vector, p_fixed, data, choice)

p = [p_fixed(1:8)' p_search']'; % here you reconstruct a proper parameter vector
sim_data= generate_ss_data(m_vector, tspan, x0, p, choice);% here you use generate_ss_data to get the appropriate data
if any(p_search<0) % tests if any parameter is negative
    cost= 2000; %typical bad costs are ~XXX
    return;
end
% Tend = tspan;
% N_comp= zeros(length(m_vector),1);
% x0= [100000 0 0 0 0];
%for i=1:length(m_vector)
%     [t,x] = ode15s(@you_odeRI, (0:Tend), x0, [], p, m_vector(i));
%     N_comp(i) = x(end,1);
% end
cost = sqrt(mean((sim_data - data).^2));
%here you compute a cost that will allow you to fit predictions to experimental data
end