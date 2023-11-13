% Mathematically define the base equations of the cycloidal gear profile
function x = base_eqn(X, dr_e, dr_p, t)                                
   Rr_e = X.r_e+dr_e; Rr_p = X.r_p+dr_p;                   
   ups = sqrt(X.E^2*(X.N+1)^2-2*X.E*Rr_e*(X.N+1)*cos(X.N*t)+Rr_e^2);  
   x = [Rr_e*(1-Rr_p./ups).*cos(t)-X.E*(1-Rr_p*(X.N+1)./ups).*cos((X.N+1)*t);
        Rr_e*(1-Rr_p./ups).*sin(t)-X.E*(1-Rr_p*(X.N+1)./ups).*sin((X.N+1)*t);];
end