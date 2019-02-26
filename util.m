% Author: Moritz Greiwe
% Description: Function that returns the utility for a given period for the
% model according to the paper of Bongini (2018) on Partial Adjustment
function util = util(z,k,kprime) % include p
global chi delta r xi tau_d tau_c A
a = length(k);
b = length(z);
k = kron(k,ones(a*b,1));
kprime = kron(ones(b,1),kron(kprime, ones(a,1)));
z = kron(ones(a*b,1),z);



%calculate profits
pi = A*z.*k.^chi;

% check for indicator with y
p =1/(1+r/(1+r))*(s*(1-delta)*k +(1-tau_c)*A*z(1)*k.^chi+tau_c*k);

% get y
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