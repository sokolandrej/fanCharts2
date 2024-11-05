function plotsofuniformity(pdens,figtitle, crit_value_vec, significance_level)

bin = 10;
m = length(pdens);

if significance_level == 1
    n_crit_value = norminv(0.01/2);
    crit_value = crit_value_vec.ks(1,1);
    lbl = '1% Critical Value';
elseif significance_level == 5
    n_crit_value = norminv(0.05/2);
    crit_value = crit_value_vec.ks(1,2);
    lbl = '5% Critical Value';
elseif significance_level == 10
    n_crit_value = norminv(0.1/2);
    crit_value = crit_value_vec.ks(1,3);
    lbl = '10% Critical Value';
else
    error('There is no critical value available for the selected significance level, choose between 1, 5 and 10');
end

% histograms
fn = figure;
hs = histc(pdens,0:1/bin:1);
bar(0:1/bin:1,bin*hs./m,'histc')
hold on
plot(0:1/bin:1,bin*(1/bin)*ones(bin+1,1), '-r','LineWidth',2)
hold on
plot(0:1/bin:1,bin*(1/bin+n_crit_value*sqrt((1/bin*(1-1/bin))/m))*ones(bin+1,1), '--r','LineWidth',2)
hold on
plot(0:1/bin:1,bin*(1/bin-n_crit_value*sqrt((1/bin*(1-1/bin))/m))*ones(bin+1,1), '--r','LineWidth',2)
xlim([0 1])
title(figtitle)
box on

% tests
rvec = [0:0.001:1];
for r = 1:size(rvec,2)
    ecdf(:,r) = mean(pdens < rvec(:,r));
end
fn = figure; 
xlabel('r') 
ylabel('\phi_P(r)') 
plot(rvec,ecdf,'LineWidth',2)
hold on
plot(rvec,rvec,'r','LineWidth',2);
hold on
plot(rvec,rvec + crit_value/sqrt(m),'r:','LineWidth',2);
hold on
plot(rvec,rvec - crit_value/sqrt(m),'r:','LineWidth',2);
hold off
xlim([0 1])
ylim([0 1])
grid on
legend('Empirical','Theoretical',lbl,'Location','NorthWest'); 
title(strcat(figtitle, ', RS (2019) Test'))
box on