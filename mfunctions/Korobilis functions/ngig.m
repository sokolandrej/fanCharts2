function [out] = ngig(lam,chi,psi)

h = lam;                   
b = sqrt( chi .* psi );

sampleSize = length(chi);

m = ( h-1+sqrt((h-1)^2 + b.^2) ) ./ b;  % Mode
 
% without shifting by m                                  
ym = (-h-1 + sqrt((h+1)^2 + b.^2))./b;  
% a = vplus/uplus  
a = exp(-0.5*h*log(m.*ym) + 0.5*log(m./ym) + b./4.*(m + 1./m - ym - 1./ym));
    
u = rand( sampleSize,1 );
v = rand( sampleSize,1 );
out = a .* (v./u);    