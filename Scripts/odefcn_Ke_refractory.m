function dydt = odefcn_Ke_refractory(t,y,beta,delta,pi,Phi,roh,k,c)

dydt = zeros(5,1);
% if t<t0
% 
%     R0 = beta*pi/(delta*c)*S0;
%     r = 1/2*(-(k+delta)+sqrt((k+delta)^2 + 4*k*delta*(R0-1)));
% 
%     ddydt(1) = -beta*y(5)*y(1); %T
%     dydt(2) = 0; %R
%     dydt(3) = beta*y(5)*y(1) - k*y(3); %E
%     dydt(4) = k*y(3) - delta*y(4); %I
%     dydt(5) = pi*y(4) - c*y(5); %V
% 
% else
    dydt(1) = -beta*y(5)*y(1) - Phi*y(4)*y(1) + roh*y(2); %T
    dydt(2) = Phi*y(4)*y(1) - roh*y(2); %R
    dydt(3) = beta*y(5)*y(1) - k*y(3); %E
    dydt(4) = k*y(3) - delta*y(4); %I
    dydt(5) = pi*y(4) - c*y(5); %V
% end

end