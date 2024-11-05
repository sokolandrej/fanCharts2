function [out] = gig(lam,chi,psi)

h = lam;                   
b = sqrt( chi * psi );
sampleSize = 1;

m = ( h-1+sqrt((h-1)^2+b^2) ) / b;  % Mode
log_1_over_pm = -(h-1)/2*log(m) + b/4*(m + (1/m));

% r = (6*m + 2*h*m - b*m^2 + b)/(4*m^2);
% s = (1 + h - b*m)/(2*m^2);
% p = (3*s - r^2)/3;
% q = (2*r^3)/27 - (r*s)/27 + b/(-4*m^2);
% eta = sqrt(-(p^3)/27);
                         
% y1  = 2*exp(log(eta)/3) * cos(acos(-q/(2*eta))/3) - r/3;
% y2  = 2*exp(log(eta)/3) * cos(acos(-q/(2*eta))/3 + 2/3*pi) - r/3;
%                          
% if (h<=1 & b<=1) | abs(q/eta)>2 | y1<0 | y2>0
    % without shifting by m                               
    ym = (-h-1 + sqrt((h+1)^2 + b^2))/b;                        
    % a = vplus/uplus
    a = exp(-0.5*h*log(m*ym) + 0.5*log(m/ym) + b/4*(m + 1/m - ym - 1.0/ym));
    
    u = rand( sampleSize );
    v = rand( sampleSize );
    out = a * (v./u);
    indxs = find( log(u) > (h-1)/2*log(out) - b/4*(out + 1./out) + log_1_over_pm );
    while ~isempty( indxs )
        indxsSize = size( indxs );
        u = rand( indxsSize );
        v = rand( indxsSize );
        outNew = a * (v./u);
        l = log(u) <= (h-1)/2*log(outNew) - b/4*(outNew + 1./outNew) + log_1_over_pm;
        out( indxs( l ) ) = outNew(l);
        indxs = indxs( ~l );
    end       
% else % if h<=1 & b<=1
%     % with shifting by m                                                  
%     vplus = exp( log_1_over_pm + log(1/y1) + (h-1)/2*log(1/y1 + m) - ...
%         b/4*(1/y1 + m + 1/(1/y1 + m)) );
%     vminus = -exp( log_1_over_pm + log(-1/y2) + (h-1)/2*log(1/y2 + m) - ...
%         b/4*(1/y2 + m + 1/(1/y2 + m)) );  
%     
%     u = rand( sampleSize );
%     v = vminus + (vplus - vminus) * rand( sampleSize );
%     z = v ./ u;
%     clear('v');
%     indxs = find( z < -m );
%     
%     while ~isempty(indxs),
%         indxsSize = size( indxs );
%         uNew = rand( indxsSize );
%         vNew = vminus + (vplus - vminus) * rand( indxsSize );
%         zNew = vNew ./ uNew;
%         l = (zNew >= -m);
%         z( indxs( l ) ) = zNew(l);
%         u( indxs( l ) ) = uNew(l);
%         indxs = indxs( ~l );
%     end
%     
%     out = z + m;   
%     indxs = find( log(u) > (log_1_over_pm + (h-1)/2*log(out) - b/4*(out + 1./out)) );                         
%    
%     while ~isempty(indxs),
%         indxsSize = size( indxs );                             
%         u = rand( indxsSize );
%         v = vminus + (vplus - vminus) * rand( indxsSize );
%         z = v ./ u;
%         clear('v');
%         indxs1 = find( z < -m );
%         while ~isempty(indxs1),
%             indxsSize1 = size( indxs1 );
%             uNew = rand( indxsSize1 );
%             vNew = vminus + (vplus - vminus) * rand( indxsSize1 );
%             zNew = vNew ./ uNew;
%             l = (zNew >= -m);
%             z( indxs1( l ) ) = zNew(l);
%             u( indxs1( l ) ) = uNew(l);
%             indxs1 = indxs1( ~l );
%         end
%         
%         outNew = z + m;
%         l = ( log(u) <= (log_1_over_pm + (h-1)/2*log(outNew) - b/4*(outNew + 1./outNew)) );
%         out( indxs(l) ) = outNew( l );
%         indxs = indxs( ~l );        
%     end        
% end %% if h<=1 & b<=1, else ...
      
out = sqrt( chi / psi ) * out;