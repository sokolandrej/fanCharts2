% OOS estimates of CQR model

%% Housekeeping and options
close all; clear;clc
addpath(genpath('mfunctions'))
quantiles = .05:.05:.95;
quantilesFC = .025:.001:.975;
quantilesToFit = quantiles([1 3 5 10 15 17 19]);
doNowcast = 0;
yTransf = 0;
horizon = 9;
nSave = 4e4;
nBurn = 1e4;
chainStep = 40;
condMode = 1;
fName = 'cqr_rep';
nVintages = 88;
lsTol = 1e-4;
outFolder = 'replication';
mkdir(outFolder);
rng(0);


%% Load and cook data
load data

% growth rates/diff yoy
data(5:end,[1:3]) = 100*log(data(5:end,[1:3])./data(1:end-4,[1:3]));
data(5:end,5) = data(5:end,5)-data(4:end-1,5);
data = data(:,1:5);

% real variables
data = [data data(:,[1 4])-data(:,2)];
data(2:end,[4 end]) =  data(2:end,[4 end])- data(1:end-1,[4 end]);
data =data(5:end,:);dates = dates(5:end);
varNames = {'nom_hp','hicp','Income','ltndiff','Unemp.', 'RHP','Int. Rate'};
data(end-horizon+1:end,2) = NaN;


%% Estimate etc
oSteps = nVintages-1;
forecasts = cell(nVintages,1);
fDates = zeros(nVintages,1);
yyFull = data(2:end,6);
dates = dates(2:end);
XXfull = [ones(length(yyFull),1) data(1:end-1,[6 3 7 5]) data(2:end,[6 3 7 5])];
XXnames = [{'const'}, strcat(varNames([6 3 7 5]),'(-2)') strcat(varNames([6 3 7 5]),'(-1)')];

% figure out conditioning horizon
Tfull = find(~isnan(sum(XXfull,2)),1,'last');
XXfull(Tfull+1:end,1:end-3) = NaN; %to avoid perfect collinearity with conditioning!
XXfull(Tfull+2:end,end-2) = NaN; %no conditioning on income save for 1q
condMask = isnan(XXfull(Tfull+1:end,:));
condHor = size(condMask,1);

if isempty(gcp('nocreate'))

    parpool
    parfevalOnAll(@warning,0,'off');

end

parfor oo = 1+(0:oSteps)

    yy = yyFull;
    XX = XXfull(1:end-oo+1,:);
    XXcondTemp = XX(end-condHor+1:end,:);
    XXcondTemp(condMask) = NaN;
    XX(end-condHor+1:end,:) = XXcondTemp;
    T0 = find(~isnan(sum(XX,2)),1,'first');
    T = find(~isnan(sum(XX,2)),1,'last');
    TC = find(sum(~isnan(XX),2),1,'last');
    yy = yy(T0:T);
    XX = XX(T0:TC,:);
    datesTemp = dates(T0:T);
    disp(['Vintage: ' num2str(datevec(datesTemp(end)))])
    lpTemp = zeros(horizon,length(quantiles));
    res = qr_cond_lp(...
        yy,XX,quantiles,horizon,XXnames,condMode,doNowcast,yTransf,nSave,nBurn,chainStep,0,[]);

    for hh = 1:horizon

        lpTemp(hh,:) = res.(['h' num2str(hh)]).lpQR';

    end

    forecasts{oo,1} = lpTemp;
    fDates(oo) = datesTemp(end);

end

forecasts = forecasts(end:-1:1);
fDates = fDates(end:-1:1);


%% Post-processing
% Fitted distributions
fittedDist = cell(nVintages,horizon);

parfor vv = 1:nVintages

    for hh = 1:horizon

        [fTemp,idx] = sort(forecasts{vv,1}(hh,:));
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

yy = data(1:end-horizon,6); %gets rid of NaNs at end
plot_lp_new(fanChart,yy(end-horizon-40:end-horizon,1),datetime(datevec(dates(end-2*horizon-41:end-horizon-1))),rgb('dark blue')); %note date shift to make axis ticks align properly in plot!
plot(datetime(datevec(dates(end-2*horizon-41:end-horizon-1))),[NaN(40,1);yy(end-horizon:end)],'r','lineWidth',2)
lgd = legend;
lgd.String{end} = 'Outturns';
xtickformat('yy-QQQ')
ylim([-8.5 8.5])
SaveFigure([outFolder '/Figure_5'],0);

% Q5 conditional density
[uncond(1), uncond(2), uncond(3), uncond(4)] = QuantilesInterpolation(quantile(data(1:end-horizon,6),quantiles),quantiles);
figure
xRange = -10:.05:10;
hh=5;
plot(xRange,dskt(xRange,fittedDist{end-horizon,hh}(1),fittedDist{end-horizon,hh}(2),fittedDist{end-horizon,hh}(3),fittedDist{end-horizon,hh}(4)),'color',rgb('dark blue'),'LineWidth',2); hold on
plot(xRange,dskt(xRange,uncond(1,1),uncond(1,2),uncond(1,3),uncond(1,4)),'--','color','k','LineWidth',2);
axis tight
yLength = get(gca,'YLim');
line([yy(end-horizon+hh),yy(end-horizon+hh)],yLength,'Color','r','LineWidth',2)
% title(datestr(dates(end-2*horizon+hh),'yy-qq'))
legend('Conditional distribution','Unconditional distribution','Outturn','Location','best')
SaveFigure([outFolder '/Figure_6a'],0);

% Q9 conditional density
figure
hh=9;
plot(xRange,dskt(xRange,fittedDist{end-horizon,hh}(1),fittedDist{end-horizon,hh}(2),fittedDist{end-horizon,hh}(3),fittedDist{end-horizon,hh}(4)),'color',rgb('dark blue'),'LineWidth',2); hold on
plot(xRange,dskt(xRange,uncond(1,1),uncond(1,2),uncond(1,3),uncond(1,4)),'--','color','k','LineWidth',2);
yLength = get(gca,'YLim');
line([yy(end-horizon+hh),yy(end-horizon+hh)],yLength,'Color','r','LineWidth',2)
% title(datestr(dates(end-2*horizon+hh),'yy-qq'))
legend('Conditional distribution','Unconditional distribution','Outturn','Location','best')
SaveFigure([outFolder '/Figure_6b'],0);

% PITs and log scores by horizon
PITs = NaN(nVintages-horizon,horizon);
qScores = NaN(nVintages-horizon,horizon,length(quantiles));
logScores = PITs;
GRScores = NaN(nVintages-horizon,horizon,5);
PITdate = PITs;
outturns = PITs;
pointFcastErr =PITs;
nBins =5;
quantilesFC = .025:.025:.975; %less fine for oos plots
quants1yAhead = NaN(nVintages-horizon,length(quantilesFC));
quants2yAhead = quants1yAhead;
quants1qAhead = quants1yAhead;

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
            fTemp = sort(forecasts{nVintages-oo,1}(hh,:));

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
shadedplot(datetime(datevec(fDates(5:end-horizon+4))),quants1yAhead(:,1)',quants1yAhead(:,end)','blue','blue')
hold on
shadedplot(datetime(datevec(fDates(5:end-horizon+4))),quants1yAhead(:,6)',quants1yAhead(:,34)','blue','blue')
hold on
plot(datetime(datevec(fDates(5:end-horizon+4))),outturns(:,4),'r','LineWidth',2)
hh = gca;
legend(hh.Children([8 4 1]),'95%','68%','Outturns','Location','best')
ylim([-15 15])
SaveFigure([outFolder '/Figure_7'],0);

figure
shadedplot(datetime(datevec(fDates(9:end-horizon+8))),quants2yAhead(:,1)',quants2yAhead(:,end)','blue','blue')
hold on
shadedplot(datetime(datevec(fDates(9:end-horizon+8))),quants2yAhead(:,6)',quants2yAhead(:,34)','blue','blue')
hold on
plot(datetime(datevec(fDates(9:end-horizon+8))),outturns(:,8),'r','LineWidth',2)
hh = gca;
legend(hh.Children([8 4 1]),'95%','68%','Outturns','Location','best')
ylim([-15 15])
SaveFigure([outFolder '/Figure_B12'],0);


%% Save outputs
save([outFolder '/' fName],'forecasts','fDates','fittedDist','PITs','logScores','qScores','GRScores','PITdate','outturns','pointFcastErr','quants1yAhead','quants2yAhead')


