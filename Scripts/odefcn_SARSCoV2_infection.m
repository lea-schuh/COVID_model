function dydt = odefcn_SARSCoV2_infection(t,y,b0,dI,pV,dV,pN,dN,pB,dB,B_thres)

    dydt = zeros(4,1);
    
    dydt(1) =  pN - dN*y(1) - b0*y(3)*y(1)*(1-y(4)) ; %S_n
    dydt(2) = b0*y(3)*y(1)*(1-y(4)) - dI*y(2); %I_n
    dydt(3) = pV*y(2) - dV*y(3) - b0*y(3)*y(1)*(1-y(4)); %V_n
    dydt(4) = pB*y(3)*(1-y(4)) - dB*(y(4)-B_thres)*y(4); %B - immune response

end