% Statistical tests of PITs and scores

%% Housekeeping and options
clear all; clc;close all
addpath(genpath('mfunctions'))
quantiles = .05:.05:.95;
outFolder = 'replication';
fNames = {[outFolder '/cqr_rep'],[outFolder '/bvar_rep']};
rng(0);

%% extract
for ff = 1:length(fNames)
    
    tempDat = load(fNames{ff},'qScores','fDates','PITs','GRScores');
    qsTemp = tempDat.qScores([1:79],:,:);
    fDates = tempDat.fDates;
    grTemp = tempDat.GRScores(1:79,:,:); 
    pTemp = tempDat.PITs([1:79],:);%46:2010Q1; 30:2006Q1 (so 26 to exclude 2008 forecast)
    avQscores{ff} = squeeze(nanmean(qsTemp,1));
    GRScores{ff} = squeeze(nanmean(grTemp,1));
    GRseries{ff} = grTemp;
    PITs{ff} = pTemp;
    
end


%% q-scores
figure

for hh = 1:size(avQscores{1},1)

    subplot(3,3,hh)
    plot (quantiles,avQscores{1}(hh,:),quantiles,avQscores{2}(hh,:),'lineWidth',2);
    ylim([0 1.2])
    if hh == 4

        ylabel('Average quantile score')

    end

    if hh == 8

            xlabel('\tau')
    end

    if hh == 1

        legend('QR','BVAR')

    end

    title(['h = ' num2str(hh)]);

end

SaveFigure([outFolder '/Figure_8'],0);


%% GR scores table
GRtab = GRScores{2}./GRScores{1};
GRpVal = GRtab;

for hh = 1:size(GRtab,1)

    for qq = 1:size(GRtab,2)

       [~,GRpVal(hh,qq)] = DM_HLNtest(GRseries{1}(:,hh,qq),GRseries{2}(:,hh,qq),hh-1,0);

    end

end

save([outFolder '/Table_2'],'GRtab','GRpVal');


%% PITs plot
figure
nBins = 20;

for hh = 1:9

    PITs_all = [PITs{1}(:,hh) PITs{2}(:,hh)];
    N = zeros(size(PITs_all,2),nBins+1);

    for nn=1:size(PITs_all,2)

        PITtemp = PITs_all(1:end,nn);
        [NTemp(nn,:),edges] = histcounts(PITtemp(~isnan(PITtemp)), nBins, 'Normalization','cdf','BinLimits',[0 1]); % CDF of PITs to compare with the Uniform CDF
        N(nn,:) = [0 NTemp(nn,:)];
        edges = edges(2:end) - (edges(2)-edges(1)); % Store edges of bars to plot a curve that fits the histogram

    end

    edges = [edges 1];
    unif_range = 0:.05:1;  % rvec for rs_test
    [~, boot_crit_value] = rs_test(PITs_all(:,1), unif_range);
    cV = boot_crit_value.ks(1,2);
    lgd_values2 = {'CQR', 'BVAR', '5% Critical Value'};
    subplot(3,3,hh)
    shadedplot([0 1],[cV/sqrt(size(PITs_all,1)-1) 1+cV/sqrt(size(PITs_all,1)-1)],[-cV/sqrt(size(PITs_all,1)-1) 1-cV/sqrt(size(PITs_all,1)-1)],'green', 'green'); 
    title(['h = ' num2str(hh)])
    hold on
    plot(edges, N(1,:),'LineWidth',2, 'Color', 'blue');
    hold on
    plot(edges, N(2,:),'LineWidth',2, 'Color', 'red');
    hold on
    plot(0:1,0:1,'Color','#77AC30','LineWidth',0.5,'LineStyle','--');
    
    if hh == 8
    
        xlabel('Quantiles');
    
    end
    
    if hh == 4
    
        ylabel('PIT')
    
    end
    
    ll = gca;
    ylim([0 1]);
    
    if hh == 9
    
        legend(ll.Children([3 2 4]), 'CQR', 'BVAR','5% C.V.','Location','best','EdgeColor','none')
    
    end

end

SaveFigure([outFolder '/Figure_B15'],0)

