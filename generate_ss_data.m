function E_ss= generate_ss_data(m_vector, tspan, x0, p, choice)
%for various values of m in m_vector computes the steady state value of the
%quantity XX and return these in the data vector
%choice is 1 for the "R" circuit, 2 for the "I" circuit, and 3 for the "RI" circuit
Tend = tspan;
E_ss = zeros(length(m_vector),1);
for m_index = 1:length(m_vector)
    m = m_vector(m_index);
    if choice==1
        [t x]= ode15s(@you_odeI, (0:Tend), x0, [], p,m);
    elseif choice==2
        [t x]= ode15s(@you_odeR, (0:Tend), x0, [], p,m);
    elseif choice==3
        [t x]= ode15s(@you_odeRI, (0:Tend), x0, [], p,m);
    else
        error('wrong choice');
    end
    %XX update the data vector
    E_ss(m_index) = x(end,1);
end
end