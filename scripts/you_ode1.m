function xdot= you_ode1(tspan,x,p)
xdot= zeros(size(x));

%N
xdot(1)= p(1)*x(1)*(1-(x(1)/p(2))) - p(3)*x(1)*x(2);
%E
xdot(2)= p(4)*x(3)*x(5) - p(5)*x(2);
%A
xdot(3)= p(6)*x(1)*x(4)  - p(7)*x(3); 
%I 
xdot(4)= p(8) - p(5)*x(4);
%R
xdot(5)= p(8)- p(5)*x(5);
end