% Author: Moritz Greiwe
% Description: Function that returns the utility for a given period for the
% model according to the paper of Bongini (2018) on Partial Adjustment
function util = util(k,p,kprime,pprime,z) % include p
global chi delta r xi tau_d tau_c A s



%calculate profits
pi = A.*z.*k.^chi;
% get y
y = pi - delta.*k -r/(1+r).*p;
% check for taxes on profits
g = tau_c.*(y>0).*y;

%check for taxes/flotation cost on cash flow
c = pi-g-p-kprime+(1-delta).*k+pprime./(1+r);
Psi = (c>0).*(1-tau_d)+(c<0).*(1+xi);


%calculate utility
util=Psi.*c+(-10000000).*(pprime > s*(1-delta)*kprime+pi-g);
end