function  out  = rho(res,tau)

res(res>0) = tau*res(res>0);
res(res<0) = (tau-1)*res(res<0);

out = res;

end

