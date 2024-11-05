function out = buildLags(X,nLags)

out = [];
if nLags>0
for ii = 1:nLags
   
%     count = count+1;
    out = [out X(nLags-ii+1:end-ii,:)];

end
end
end