function [ se ] = deltalw_ce( g,rp )
% Delta method: Ledoit Wolf
% input: g = gamma/2 
%        rp = Tx2 matrix of out-of-sample portfolio returns 
   % the CE difference will be computed as rp(:,1)-rp(:,2)
T = length(rp);
ri = rp(:,1);
rn = rp(:,2);
ri2 = ri.^2;
rn2 = rn.^2;

y = [ri-mean(ri),rn-mean(rn), ri2-mean(ri2),rn2-mean(rn2)];

psi = (1/T)*(y'*y);
der = [1+2*g*mean(ri),-1-2*g*mean(rn),-g,g]';

se = sqrt((der'*psi*der)/T);


end

