% Author: Moritz Greiwe
% Description: Function that returns the utility for a given period for the
% model according to the paper of Bongini (2018) on Partial Adjustment
function util = util(z,k,p,kprime, pprime)
global chi delta r xi tau_d tau_c A
k = kron(k, ones(length(k),1));
kprime = kron(ones(length(k),1), k);
z = kron(z, ones(length(k),1));

%calculate profits
pi = A*z.*k.^chi;
y = pi - delta*k -r/(1+r)*p;

% check for taxes on profits
indicator = y>0
g = tau_c*indicator.*y;

%check for taxes/flotation cost on cash flow
c = pi-g-p-kprime+(1-delta).*k+pprime./(1+r);
Psi1 = (c>0)*(1-tau_d);
Psi2 = (c<0)*(1+xi);
Psi = Psi1+Psi2

%calculate utility
util=Psi.*c
end