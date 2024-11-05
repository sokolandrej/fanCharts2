% OOS estimates of BVAR model

%% Housekeeping and options
close all; clear;clc
addpath(genpath('mfunctions'))
quantiles = .05:.05:.95;
quantilesFC = .025:.001:.975;
doNowcast = 0;
yTransf = 0;
horizon = 9;
nDraws = 1000;
lags = 5;
fName = 'bvar_rep';
conditioning = 1;
doParallel = 1;
nVintages = 88;
lsTol = 1e-4;
outFolder = 'replication';
mkdir(outFolder);


%% Load and cook data
load data
data(:,[1:3]) = log(data(:,[1:3]));
data = data(:,1:5);

% real variables
data = [data data(:,1)-data(:,2)];
data = [data(2:end,:) (data(2:end,4) - 400*(data(2:end,2)-data(1:end-1,2)))];
dates = dates(2:end);
varNames = {'nom_hp','hicp','pyr','ltn','ur', 'real_hp','real_lt'};
data(end-horizon+1:end,2) = NaN;


%% Estimate etc
oSteps = nVintages-1;
forecasts = cell(nVintages,1);
XXfull = data(:,[6 3 7 5]);
XXnames = varNames([6 3 7 5]);

% figure out conditioning horizon
Tfull = find(~isnan(sum(XXfull,2)),1,'last');
XXfull(Tfull+2:end,2) = NaN; %no conditioning on income save for 1q
condMask = isnan(XXfull(Tfull+1:end,:));
condHor = size(condMask,1);

if isempty(gcp('nocreate'))

    parpool
    parfevalOnAll(@warning,0,'off');

end

parfor oo = 1+(0:oSteps)

    XX = XXfull(1:end-oo+1,:);
    XXcondTemp = XX(end-condHor+1:end,:);
    XXcondTemp(condMask) = NaN;
    XX(end-condHor+1:end,:) = XXcondTemp;
    T0 = find(~isnan(sum(XX,2)),1,'first');
    T = find(~isnan(sum(XX,2)),1,'last');
    TC = find(sum(~isnan(XX),2),1,'last');
    XX = XX(T0:TC,:);

    if ~conditioning

        XX(TC-horizon+1:end,:) = NaN;

    end

    datesTemp = dates(T0:T);
    disp(['Vintage: ' num2str(datevec(datesTemp(end)))])
    res = qbvar_alt(XX,lags,[],.2,1,1,2,0,nDraws,nDraws);
    yDraws = squeeze(res.yMCMC(:,1,:));
    yDraws(5:end,:) = 100*(yDraws(5:end,:)-yDraws(1:end-4,:));
    yDraws = yDraws(end-horizon+1:end,:);
    forecasts{oo} = yDraws;
    fDates(oo) = datesTemp(end);

end

forecasts = forecasts(end:-1:1);
fDates = fDates(end:-1:1);


%% Post-processing
% Fitted distributions
fittedDist = cell(nVintages,horizon);

parfor vv = 1:nVintages

    for hh = 1:horizon

        tempData = forecasts{vv}(hh,:);
        fTemp = quantile(tempData,quantiles);
        params = zeros(4,1);
        [params(1), params(2), params(3), params(4)] = QuantilesInterpolation(fTemp,quantiles);
        fittedDist{vv,hh} = params;

    end

end

% Latest fanchart with full set of outturns
fanChart = zeros(horizon,length(quantilesFC));

for hh = 1:horizon

    if ~isempty(fittedDist{end-horizon,hh})

        fanChart(hh,:) = qskt(quantilesFC,fittedDist{end-horizon,hh}(1),fittedDist{end-horizon,hh}(2),fittedDist{end-horizon,hh}(3),fittedDist{end-horizon,hh}(4));

    else

        fanChart(hh,:) = zeros(1,length(quantilesFC));

    end

end

yy = 100*(data(5:end,6)-data(1:end-4,6));
yy = yy(~isnan(yy));
plot_lp_new(fanChart,yy(end-horizon-40:end-horizon,1),datetime(datevec(dates(end-2*horizon-41:end-horizon-1))),rgb('orange'));
plot(datetime(datevec(dates(end-2*horizon-41:end-horizon-1))),[NaN(40,1);yy(end-horizon:end)],'r','lineWidth',2)
ylim([-8.5 8.5])
lgd = legend;
lgd.String{end} = 'Outturns';
xtickformat('yy-QQQ')
SaveFigure([outFolder '/Figure_B13'],0);

% PITs and log scores by horizon
PITs = NaN(nVintages-horizon,horizon);
quantiles = .05:.05:.95;
qScores = NaN(nVintages-horizon,horizon,length(quantiles));
GRScores = NaN(nVintages-horizon,horizon,5);
logScores = PITs;
PITdate = PITs;
outturns = PITs;
pointFcastErr =PITs;
nBins =5;
quantilesFC = .025:.025:.975; %less fine for oos plots
quants1yAhead = NaN(nVintages-horizon,length(quantilesFC));
quants2yAhead = quants1yAhead;


for hh = 1:horizon

    for oo = nVintages-1:-1:horizon

        if ~isempty(fittedDist{nVintages-oo,hh})

            GRScores(nVintages-oo,hh,:) = skewt_GRscores(yy(end-oo+hh),[fittedDist{nVintages-oo,hh}(1),fittedDist{nVintages-oo,hh}(2),fittedDist{nVintages-oo,hh}(3),fittedDist{nVintages-oo,hh}(4)]);
            PITs(nVintages-oo,hh) = pskt(yy(end-oo+hh),fittedDist{nVintages-oo,hh}(1),fittedDist{nVintages-oo,hh}(2),fittedDist{nVintages-oo,hh}(3),fittedDist{nVintages-oo,hh}(4));
            tempLogScore = dskt(yy(end-oo+hh),fittedDist{nVintages-oo,hh}(1),fittedDist{nVintages-oo,hh}(2),fittedDist{nVintages-oo,hh}(3),fittedDist{nVintages-oo,hh}(4));
            tempLogScore(tempLogScore<lsTol) = lsTol;
            logScores(nVintages-oo,hh) = log(tempLogScore);
            outturns(nVintages-oo,hh) = yy(end-oo+hh);
            pointFcastErr(nVintages-oo,hh) =  qskt(.5,fittedDist{nVintages-oo,hh}(1),fittedDist{nVintages-oo,hh}(2),fittedDist{nVintages-oo,hh}(3),fittedDist{nVintages-oo,hh}(4),1e-6) - yy(end-oo+hh);
            fTemp = qskt(quantiles,fittedDist{nVintages-oo,hh}(1),fittedDist{nVintages-oo,hh}(2),fittedDist{nVintages-oo,hh}(3),fittedDist{nVintages-oo,hh}(4),1e-6);

            for qq = 1:length(quantiles)

                qScores(nVintages-oo,hh,qq) = rho(outturns(nVintages-oo,hh)-fTemp(qq),quantiles(qq));

            end


            if hh == 4

                quants1yAhead(nVintages-oo,:) = qskt(quantilesFC,fittedDist{nVintages-oo,4}(1),fittedDist{nVintages-oo,4}(2),fittedDist{nVintages-oo,4}(3),fittedDist{nVintages-oo,4}(4),1e-6);

            elseif hh == 8

                quants2yAhead(nVintages-oo,:) = qskt(quantilesFC,fittedDist{nVintages-oo,8}(1),fittedDist{nVintages-oo,8}(2),fittedDist{nVintages-oo,8}(3),fittedDist{nVintages-oo,8}(4),1e-6);

            elseif hh == 1

                quants1qAhead(nVintages-oo,:) = qskt(quantilesFC,fittedDist{nVintages-oo,1}(1),fittedDist{nVintages-oo,1}(2),fittedDist{nVintages-oo,1}(3),fittedDist{nVintages-oo,1}(4),1e-6);


            end

        end

        PITdate(nVintages-oo,hh) = datenum(dates(end-horizon-oo+hh));

    end

end

figure
shadedplot(datetime(datevec(fDates(5:end-horizon+4))),quants1yAhead(:,1)',quants1yAhead(:,end)','red','red')
hold on
shadedplot(datetime(datevec(fDates(5:end-horizon+4))),quants1yAhead(:,6)',quants1yAhead(:,34)','red','red')
hold on
plot(datetime(datevec(fDates(5:end-horizon+4))),outturns(:,4),'b','LineWidth',2)
ff = gca;
legend(ff.Children([8 4 1]),'95%','68%','Outturns')
SaveFigure([outFolder '/Figure_B14a'],0);

figure
shadedplot(datetime(datevec(fDates(9:end-horizon+8))),quants2yAhead(:,1)',quants2yAhead(:,end)','red','red')
hold on
shadedplot(datetime(datevec(fDates(9:end-horizon+8))),quants2yAhead(:,6)',quants2yAhead(:,34)','red','red')
hold on
plot(datetime(datevec(fDates(9:end-horizon+8))),outturns(:,8),'b','LineWidth',2)
ff = gca;
legend(ff.Children([8 4 1]),'95%','68%','Outturns')
SaveFigure([outFolder '/Figure_B14b'],0);

save([outFolder '/' fName],'forecasts','fDates','fittedDist','PITs','GRScores','logScores','qScores','PITdate','outturns','pointFcastErr','quants1yAhead','quants2yAhead')
