% scenarios for paper

%% Housekeeping and options
close all; clear;clc
addpath(genpath('mfunctions'))
quantiles = .05:.05:.95;
quantilesFC = .025:.001:.975;%[.025 .16 .4999 .5 .5001 .84 .975];
quantilesToFit = quantiles([1 3 5 10 15 17 19]);%[.05 .15 .5 .85]; %weird rounding problem...
doNowcast = 0;
yTransf = 0; %need to change units of observed data if doing this
horizon = 9;
nSave = 4e4;
nBurn = 1e4;
chainStep = 40;
condMode = 2; 
fName = 'scenarios_rep';
lsTol = 1e-4;
outFolder = 'replication';
mkdir(outFolder);
rng(0);

%% Load and cook data
load data
% growth rates/diff yoy
data(5:end,[1:3]) = 100*log(data(5:end,[1:3])./data(1:end-4,[1:3]));
data = data(:,1:5);

% real variables
data = [data data(:,[1 4])-data(:,2)];
dForCond = data;
data(5:end,5) = data(5:end,5)-data(4:end-1,5);
dForCond = dForCond(5:end,:);
data(2:end,[4 end]) =  data(2:end,[4 end])- data(1:end-1,[4 end]);
data =data(5:end,:);dates = dates(5:end);
varNames = {'nom_hp','hicp','Income','ltndiff','Unemp.', 'RHP','Int. Rate'};
data(end-horizon+2:end,2) = NaN;


%% Estimate etc
yy = data(2:end,6);
dates = dates(2:end);
XX = [ones(length(yy),1) data(1:end-1,[6 3 7 5]) data(2:end,[6 3 7 5])];
XXnames = [{'Const'}, strcat(varNames([6 3 7 5]),'(-1)') varNames([6 3 7 5])];

% figure out conditioning horizon
Tfull = find(~isnan(sum(XX,2)),1,'last'); 
XX(Tfull+1:end,1:end-3) = NaN;
XX(Tfull+2:end,end-2) = NaN;
condMask = isnan(XX(Tfull+1:end,:));
condHor = size(condMask,1);
condVar = ~isnan(XX(end,:));

%first version: flat (no growth) conditioning paths
XX1 = XX;
XX1(Tfull+2:end,condVar) = repmat(0*XX1(Tfull+1,condVar),horizon-1,1); %because house prices lag everything else...
res1 = qr_cond_lp(...
    yy(1:Tfull),XX1,quantiles,horizon,XXnames,condMode,doNowcast,yTransf,nSave,nBurn,chainStep,1,[]);

for hh = 1:horizon
    
    lp1(hh,:) = res1.(['h' num2str(hh)]).lpQR';
    
end

%second version: recovery: unemp 1% lower in 2 years, LTR 1% lower
XX2 = XX;
XX2(Tfull+2:end,condVar) = repmat([-.125 -.125],horizon-1,1);

res2 = qr_cond_lp(...
    yy(1:Tfull),XX2,quantiles,horizon,XXnames,condMode,doNowcast,yTransf,nSave,nBurn,chainStep,1,[]);

for hh = 1:horizon
    
    lp2(hh,:) = res2.(['h' num2str(hh)]).lpQR';
    
end

%third version: use info up to T only
XX3 = XX;
XX3(Tfull+2:end,condVar) = NaN; %because house prices lag everything else...


res3 = qr_cond_lp(...
    yy(1:Tfull),XX3,quantiles,horizon,XXnames,condMode,doNowcast,yTransf,nSave,nBurn,chainStep,1,[]); 

for hh = 1:horizon
    
    lp3(hh,:) = res3.(['h' num2str(hh)]).lpQR';
    
end

%fourth version: stress: unemp 2% higher in 2 years, LTR 1% higher
XX4 = XX;
XX4(Tfull+2:end,condVar) = repmat([.125 .25],horizon-1,1); %because house prices lag everything else...

res4 = qr_cond_lp(...
    yy(1:Tfull),XX4,quantiles,horizon,XXnames,condMode,doNowcast,yTransf,nSave,nBurn,chainStep,1,[]);

for hh = 1:horizon
    
    lp4(hh,:) = res4.(['h' num2str(hh)]).lpQR';
    
end


%% Plots etc
% Coefficient plots
fnCount = 1;

for hh = [1 4 8]

    draws = res1.(['h' num2str(hh)]).bQRdraws;
    labels = res1.(['h' num2str(hh)]).XXlabels;
    labels = strrep(labels,'_',' ');
    bb =mean(draws,3);
    bL = quantile(draws,.16,3);
    bH = quantile(draws,.84,3);
    fnCount = fnCount + 1;
    plot_qr_coefficients(bb,bL,bH,[],[],[],quantiles,68,labels,'blue',[outFolder '/Figure_' num2str(fnCount)]);

end

% Conditioning variables
[~,cIdx] = find(condVar);
uStart = dForCond(end-horizon-41,5);
rStart = dForCond(end-horizon-41,end);
levStart = [rStart uStart];
suffix = {'a','b'};

for cc = 1:sum(condVar)
    
    figure
    plot(datetime(datevec(dates(end-horizon-41:end-1))),levStart(cc)+cumsum(XX1(end-horizon-40:end,cIdx(cc))),'color',rgb('blue'),'lineWidth',2);
    hold on;
    plot(datetime(datevec(dates(end-horizon-41:end-1))),levStart(cc)+cumsum(XX2(end-horizon-40:end,cIdx(cc))),'color',rgb('green'),'lineWidth',2);
    plot(datetime(datevec(dates(end-horizon-41:end-1))),levStart(cc)+cumsum(XX4(end-horizon-40:end,cIdx(cc))),'color',rgb('red'),'lineWidth',2);
    plot(datetime(datevec(dates(end-horizon-41:end-1))),levStart(cc)+cumsum(XX3(end-horizon-40:end,cIdx(cc))),'color',rgb('black'),'lineWidth',2);   
    legend('No change','Recovery','Stress','Available data','','Location','best');
    SaveFigure([outFolder '/Figure_9' suffix{cc}],0);
    
end

% Fitted distributions
fittedDist = cell(4,horizon);

for hh = 1:horizon
    
    [fTemp1,idx1] = sort(lp1(hh,:));
    [fTemp2,idx2] = sort(lp2(hh,:));
    [fTemp3,idx3] = sort(lp3(hh,:));
    [fTemp4,idx4] = sort(lp4(hh,:));
    params1 = zeros(4,1);
    [params1(1), params1(2), params1(3), params1(4)] = QuantilesInterpolation(fTemp1,quantiles);
    fittedDist{1,hh} = params1;
    params2 = zeros(4,1);
    [params2(1), params2(2), params2(3), params2(4)] = QuantilesInterpolation(fTemp2,quantiles);
    fittedDist{2,hh} = params2;
    params3 = zeros(4,1);
    [params3(1), params3(2), params3(3), params3(4)] = QuantilesInterpolation(fTemp3,quantiles);
    fittedDist{3,hh} = params3;
    params4 = zeros(4,1);
    [params4(1), params4(2), params4(3), params4(4)] = QuantilesInterpolation(fTemp4,quantiles);
    fittedDist{4,hh} = params4;
   
end

% Fancharts
fanChart1 = zeros(horizon,length(quantilesFC));

parfor hh = 1:horizon
    
    if ~isempty(fittedDist{1,hh})
        
        fanChart1(hh,:) = qskt(quantilesFC,fittedDist{1,hh}(1),fittedDist{1,hh}(2),fittedDist{1,hh}(3),fittedDist{1,hh}(4));
        
    else
        
        fanChart1(hh,:) = zeros(1,length(quantilesFC));
        
    end
end

fanChart2 = zeros(horizon,length(quantilesFC));

parfor hh = 1:horizon
    
    if ~isempty(fittedDist{2,hh})
        
        fanChart2(hh,:) = qskt(quantilesFC,fittedDist{2,hh}(1),fittedDist{2,hh}(2),fittedDist{2,hh}(3),fittedDist{2,hh}(4));
        
    else
        
        fanChart2(hh,:) = zeros(1,length(quantilesFC));
        
    end

end

fanChart3 = zeros(horizon,length(quantilesFC));

parfor hh = 1:horizon
    
    if ~isempty(fittedDist{3,hh})
        
        fanChart3(hh,:) = qskt(quantilesFC,fittedDist{3,hh}(1),fittedDist{3,hh}(2),fittedDist{3,hh}(3),fittedDist{3,hh}(4));
        
    else
        
        fanChart3(hh,:) = zeros(1,length(quantilesFC));
        
    end

end

fanChart4 = zeros(horizon,length(quantilesFC));

parfor hh = 1:horizon
    
    if ~isempty(fittedDist{4,hh})
        
        fanChart4(hh,:) = qskt(quantilesFC,fittedDist{4,hh}(1),fittedDist{4,hh}(2),fittedDist{4,hh}(3),fittedDist{4,hh}(4));
        
    else
        
        fanChart4(hh,:) = zeros(1,length(quantilesFC));
        
    end

end

yy = data(1:end-horizon,6); %gets rid of NaNs at end
plot_lp_new(fanChart1,yy(end-40:end,1),datetime(datevec(dates(end-horizon-41:end-1))),rgb('blue')); %note stupid date shift to make axis ticks align properly!
xtickformat('yy-QQQ')
ylim([-8.5 8.5])
SaveFigure([outFolder '/Figure_10b'],0);
plot_lp_new(fanChart2,yy(end-40:end,1),datetime(datevec(dates(end-horizon-41:end-1))),rgb('green')); %note stupid date shift to make axis ticks align properly!
xtickformat('yy-QQQ')
ylim([-8.5 8.5])
SaveFigure([outFolder '/Figure_10c'],0);
plot_lp_new(fanChart3,yy(end-40:end,1),datetime(datevec(dates(end-horizon-41:end-1))),rgb('black')); %note stupid date shift to make axis ticks align properly!
xtickformat('yy-QQQ')
ylim([-8.5 8.5])
SaveFigure([outFolder '/Figure_10a'],0);
SaveFigure([outFolder '/Figure_1a'],0);
plot_lp_new(fanChart4,yy(end-40:end,1),datetime(datevec(dates(end-horizon-41:end-1))),rgb('red')); %note stupid date shift to make axis ticks align properly!
xtickformat('yy-QQQ')
ylim([-8.5 8.5])
SaveFigure([outFolder '/Figure_10d'],0);
SaveFigure([outFolder '/Figure_1b'],0);

% Q4 and Q8 conditional densities
[uncond(1), uncond(2), uncond(3), uncond(4)] = QuantilesInterpolation(quantile(data(1:end-horizon,6),quantiles),quantiles);
xRange = -10:.05:10;
counter = 0;

for hh = [5 9]
    
    counter = counter+1;
    figure
    plot(xRange,dskt(xRange,fittedDist{1,hh}(1),fittedDist{1,hh}(2),fittedDist{1,hh}(3),fittedDist{1,hh}(4)),'color',rgb('blue'),'LineWidth',2); hold on
    plot(xRange,dskt(xRange,fittedDist{2,hh}(1),fittedDist{2,hh}(2),fittedDist{2,hh}(3),fittedDist{2,hh}(4)),'color',rgb('green'),'LineWidth',2);
    plot(xRange,dskt(xRange,fittedDist{3,hh}(1),fittedDist{3,hh}(2),fittedDist{3,hh}(3),fittedDist{3,hh}(4)),'color',rgb('black'),'LineWidth',2);
    plot(xRange,dskt(xRange,fittedDist{4,hh}(1),fittedDist{4,hh}(2),fittedDist{4,hh}(3),fittedDist{4,hh}(4)),'color',rgb('red'),'LineWidth',2);
    plot(xRange,dskt(xRange,uncond(1,1),uncond(1,2),uncond(1,3),uncond(1,4)),'--','color','k','LineWidth',2);
    legend('No change','Recovery','Available data','Stress','Unconditional distribution','Location','best')
    SaveFigure([outFolder '/Figure_11' suffix{counter}],0);
    
end

