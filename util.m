% Author: Moritz Greiwe
% Description: Function that returns the utility for a given period for the
% model according to the paper of Bongini (2018) on Partial Adjustment
function util = util(z,k,p,kprime, pprime)
global chi delta r xi tau_d tau_c A
pi = A*z*k^chi;
y = pi - delta*k -r/(1+r)*p;
if y>0
    g = tau_c*y;
else
    g=0;
end
c = pi-g-p-kprime+(1-delta)*k+pprime/(1+r);
if c>0
    Psi = 1-tau_d
else
    Psi = 1+ xi
end
util=Psi*c
end