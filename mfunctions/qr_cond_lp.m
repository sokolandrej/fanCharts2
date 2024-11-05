function fOutput = qr_cond_lp(...
    y,X,tau,hz,labels,condMode,doNowcast,yTransf,nSave,nBurn,chainStep,doParallel,bQRoptions)
%% qr_cond_lp_bma
% COMPUTES QUANTILE LOCAL PROJECTION (DIRECT FORECAST) USING BAYESIAN QUANTILE REGRESSION
% WITH A SPIKE-AND-SLAB PRIOR ON THE COEFFICIENTS, WITH OPTIONAL CONDITIONAL INFORMATION.
% condMode options:
%   1: uses all conditioning information at horizon h
%   2: uses all conditioning information up to horizon h
%   3: uses all conditioning information up to maximum conditioning horizon


if length(hz) == 1
    
    horizon = hz;
    counter = 0:hz;
    
else
    
    horizon = length(hz);
    
    if hz(1) == 0
        
        counter = hz;
        horizon = horizon-1;
        
    else
        
        counter = [0 hz];
        
    end
    
end

% parfor admin
hhAll = 1-doNowcast:horizon;
hhL = length(hhAll);
bQRdraws = cell(hhL,1);
lpQR = bQRdraws;
lpQRdraws = lpQR;
forProjectionStore = lpQR;
XXstore = lpQR;
yyStore = lpQR;
bQR = lpQR;
XXlabels = bQR;
qXingStore = bQR;

if condMode == 0
    
    X = X(1:length(y),:);
    condHor = 0;
    condIdx = [];
    labelsCond = [];
    condHorStartIdx = [];
    
end

if size(X,1)>length(y)
    
    condHor = sum(~isnan(X(length(y)+1:end,:)),1);
    condIdx = find(condHor>0);
    
    if length(condHor) >1
        
        condHorStartIdx = condHor;
        condHorStartIdx(condIdx(1)) = 1;
        
        for cc = 2:length(condIdx)
            
            condHorStartIdx(condIdx(cc)) =  condHorStartIdx(condIdx(cc-1))+condHor(condIdx(cc-1));
            
        end
        
    else
        
        condHorStartIdx = 1;
        
    end
    
    labelsCond = labels(condIdx);
    
end

if doParallel
    
    parfor hhh = 1:hhL
        
        hh = hhAll(hhh);
        
        disp(['Working on projection horizon ' num2str(counter(hh+1))])
        
        if yTransf == 1 && counter(hh+1)>0
            
            Xma = movsum(X,[counter(hh+1)-1 0],1)/counter(hh+1);
            
        else
            
            Xma = X;
            
        end
        
        if size(X,1)>length(y)
            
            if condMode == 3
                
                yEstim = y(1:end-max(condHor)+min(counter(hh+1),max(condHor)));
                Xestim = Xma(1:length(yEstim),:);
                Xcond = zeros(size(Xestim,1),sum(condHor));
                XcondProj = zeros(1,size(Xcond,2));
                labelsTemp = {};
                
                for cc = 1:length(condIdx)
                    
                    labelsTemp = [labelsTemp strcat(labelsCond(cc),strcat(strcat('(+',compose('%d',condHor(condIdx(cc)):-1:1)),')'))];
                    lagTemp = buildLags(Xma(1:length(yEstim)+condHor(condIdx(cc)),condIdx(cc)),condHor(condIdx(cc))-1);
                    XcondTemp = [Xma(condHor(condIdx(cc))+1:length(yEstim)+condHor(condIdx(cc)),condIdx(cc)) lagTemp(2:end,:)];
                    Xcond(:,condHorStartIdx(condIdx(cc)):sum(condHor(1:condIdx(cc)))) = XcondTemp(end-length(yEstim)+1:end,:);
                    
                    lagTempProj = buildLags(Xma(length(y)+1:length(y)+condHor(condIdx(cc)),condIdx(cc)),condHor(condIdx(cc))-1);
                    
                    if isempty(lagTempProj)
                        
                        XcondTempProj = Xma(length(y)+condHor(condIdx(cc)),condIdx(cc));
                        
                    else
                        
                        XcondTempProj = [Xma(length(y)+condHor(condIdx(cc)),condIdx(cc)) lagTempProj(end,:)];
                        
                    end
                    
                    XcondProj(:,condHorStartIdx(condIdx(cc)):sum(condHor(1:condIdx(cc)))) = XcondTempProj;
                    
                end
                
            else
                
                condHorTemp = min(condHor,counter(hh+1));
                
                if length(condHorTemp) >1
                    
                    condHorStartIdxTemp = condHorTemp;
                    condHorStartIdxTemp(condIdx(1)) = 1;
                    
                    for cc = 2:length(condIdx)
                        
                        condHorStartIdxTemp(condIdx(cc)) =  condHorStartIdxTemp(condIdx(cc-1))+condHorTemp(condIdx(cc-1));
                        
                    end
                    
                else
                    
                    condHorStartIdxTemp = 1;
                    
                end
                
                yEstim = y;
                Xestim = Xma(1:length(yEstim),:);
                Xcond = zeros(size(Xestim,1),sum(condHorTemp));
                XcondProj = zeros(1,size(Xcond,2));
                labelsTemp = {};
                
                for cc = 1:length(condIdx)
                    
                    labelsTemp = [labelsTemp strcat(labelsCond(cc),strcat(strcat('(+',compose('%d',condHorTemp(condIdx(cc)):-1:1)),')'))];
                    lagTemp = buildLags(Xma(1:length(yEstim)+condHorTemp(condIdx(cc)),condIdx(cc)),condHorTemp(condIdx(cc))-1);
                    XcondTemp = [Xma(condHorTemp(condIdx(cc))+1:length(yEstim)+condHorTemp(condIdx(cc)),condIdx(cc)) lagTemp(2:end,:)];
                    Xcond(:,condHorStartIdxTemp(condIdx(cc)):sum(condHorTemp(1:condIdx(cc)))) = XcondTemp(end-length(yEstim)+1:end,:);
                    
                    lagTempProj = buildLags(Xma(length(y)+1:length(y)+condHorTemp(condIdx(cc)),condIdx(cc)),condHorTemp(condIdx(cc))-1);
                    
                    if isempty(lagTempProj)
                        
                        XcondTempProj = Xma(length(y)+condHorTemp(condIdx(cc)),condIdx(cc));
                        
                    else
                        
                        
                        XcondTempProj = [Xma(length(y)+condHorTemp(condIdx(cc)),condIdx(cc)) lagTempProj(end,:)];
                        
                    end
                    
                    XcondProj(:,condHorStartIdxTemp(condIdx(cc)):sum(condHorTemp(1:condIdx(cc)))) = XcondTempProj;
                    
                end
                
                if condMode == 1
                    
                    Xcond = Xcond(:,condHorStartIdxTemp(condIdx));
                    XcondProj = XcondProj(:,condHorStartIdxTemp(condIdx));
                    labelsTemp = labelsTemp(condHorStartIdxTemp(condIdx));
                    
                end
                
            end
            
        else
            
            XcondProj = [];
            Xestim = Xma(1:length(y),:);
            yEstim = y;
            labelsTemp = {};
            
        end
        
        if yTransf == 1 && counter(hh+1)>0
            
            yy = movsum(yEstim,[0 counter(hh+1)])-yEstim;
            yy = yy(1:end-counter(hh+1))/counter(hh+1);
            
            
        else
            
            yy = yEstim(1+counter(hh+1):end);
            
        end
        
        if isempty(XcondProj)
            
            XX = Xestim(1:end-counter(hh+1),:);
            
        else
            
            XX = [Xestim(1:end-counter(hh+1),:) Xcond(1:end-counter(hh+1),:)];
            
        end
        
        bQRdraws{hhh} = qreg_bayesian( yy, XX,tau,nSave,nBurn,chainStep,bQRoptions);
        %         bQRtemp = mean(bQRdraws{hhh},3);
        bQRtemp = bQRdraws{hhh};
        %         bQR{hhh} = bQRtemp;
        bQR{hhh} = mean(bQRtemp,3);
        lpQRtemp = zeros(length(tau),nSave/chainStep);
        qXing = zeros(nSave/chainStep,1);
        forProjection = [X(length(y),:) XcondProj];
        
        for tt = 1:length(tau)

            lpQRtemp(tt,:) = ...
                forProjection*squeeze(bQRtemp(:,tt,:));

            for dd = 1:nSave/chainStep

                [~,idx] = sort(lpQRtemp(:,dd));

                if ~isequal(idx',1:length(tau))

                    qXing(dd) = 1;

                end

            end

        end
        
        %         lpQR{hhh} = mean(lpQRtemp(:,logical(1-qXing)),2);
        lpQR{hhh} = mean(lpQRtemp,2);
        lpQRdraws{hhh} = lpQRtemp;
        forProjectionStore{hhh} = forProjection;
        XXstore{hhh} = XX;
        yyStore{hhh} = yy;
        XXlabels{hhh} = [labels labelsTemp];
        qXingStore{hhh} =qXing;
        
    end
    
else %below not updated for qXing...
    
    for hhh = 1:hhL
        
        hh = hhAll(hhh);
        
%         disp(['Working on projection horizon ' num2str(counter(hh+1))])
        
        if size(X,1)>length(y)
            
            if condMode == 3
                
                yEstim = y(1:end-max(condHor)+min(counter(hh+1),max(condHor)));
                Xestim = X(1:length(yEstim),:);
                Xcond = zeros(size(Xestim,1),sum(condHor));
                XcondProj = zeros(1,size(Xcond,2));
                labelsTemp = {};
                
                for cc = 1:length(condIdx)
                    
                    labelsTemp = [labelsTemp strcat(labelsCond(cc),strcat(strcat('(+',compose('%d',condHor(condIdx(cc)):-1:1)),')'))];
                    lagTemp = buildLags(X(1:length(yEstim)+condHor(condIdx(cc)),condIdx(cc)),condHor(condIdx(cc))-1);
                    XcondTemp = [X(condHor(condIdx(cc))+1:length(yEstim)+condHor(condIdx(cc)),condIdx(cc)) lagTemp(2:end,:)];
                    Xcond(:,condHorStartIdx(condIdx(cc)):sum(condHor(1:condIdx(cc)))) = XcondTemp(end-length(yEstim)+1:end,:);
                    
                    lagTempProj = buildLags(X(length(y)+1:length(y)+condHor(condIdx(cc)),condIdx(cc)),condHor(condIdx(cc))-1);
                    
                    if isempty(lagTempProj)
                        
                        XcondTempProj = X(length(y)+condHor(condIdx(cc)),condIdx(cc));
                        
                    else
                        
                        XcondTempProj = [X(length(y)+condHor(condIdx(cc)),condIdx(cc)) lagTempProj(end,:)];
                        
                    end
                    XcondProj(:,condHorStartIdx(condIdx(cc)):sum(condHor(1:condIdx(cc)))) = XcondTempProj;
                    
                end
                
            else
                
                condHorTemp = min(condHor,counter(hh+1));
                
                if length(condHorTemp) >1
                    
                    
                    condHorStartIdxTemp = condHorTemp;
                    condHorStartIdxTemp(condIdx(1)) = 1;
                    
                    for cc = 2:length(condIdx)
                        
                        condHorStartIdxTemp(condIdx(cc)) =  condHorStartIdxTemp(condIdx(cc-1))+condHorTemp(condIdx(cc-1));
                        
                    end
                    
                else
                    
                    condHorStartIdxTemp = 1;
                    
                end
                
                yEstim = y;
                Xestim = X(1:length(yEstim),:);
                Xcond = zeros(size(Xestim,1),sum(condHorTemp));
                XcondProj = zeros(1,size(Xcond,2));
                labelsTemp = {};
                
                for cc = 1:length(condIdx)
                    
                    labelsTemp = [labelsTemp strcat(labelsCond(cc),strcat(strcat('(+',compose('%d',condHorTemp(condIdx(cc)):-1:1)),')'))];
                    lagTemp = buildLags(X(1:length(yEstim)+condHorTemp(condIdx(cc)),condIdx(cc)),condHorTemp(condIdx(cc))-1);
                    XcondTemp = [X(condHorTemp(condIdx(cc))+1:length(yEstim)+condHorTemp(condIdx(cc)),condIdx(cc)) lagTemp(2:end,:)];
                    Xcond(:,condHorStartIdxTemp(condIdx(cc)):sum(condHorTemp(1:condIdx(cc)))) = XcondTemp(end-length(yEstim)+1:end,:);
                    
                    lagTempProj = buildLags(X(length(y)+1:length(y)+condHorTemp(condIdx(cc)),condIdx(cc)),condHorTemp(condIdx(cc))-1);
                    
                    if isempty(lagTempProj)
                        
                        XcondTempProj = X(length(y)+condHorTemp(condIdx(cc)),condIdx(cc));
                        
                    else
                        
                        
                        XcondTempProj = [X(length(y)+condHorTemp(condIdx(cc)),condIdx(cc)) lagTempProj(end,:)];
                        
                    end
                    
                    XcondProj(:,condHorStartIdxTemp(condIdx(cc)):sum(condHorTemp(1:condIdx(cc)))) = XcondTempProj;
                    
                end
                
                if condMode == 1
                    
                    Xcond = Xcond(:,condHorStartIdxTemp(condIdx));
                    XcondProj = XcondProj(:,condHorStartIdxTemp(condIdx));
                    labelsTemp = labelsTemp(condHorStartIdxTemp(condIdx));
                    
                end
                
            end
            
        else
            
            XcondProj = [];
            Xestim = X(1:length(y),:);
            yEstim = y;
            labelsTemp = {};
            
        end
        
        if yTransf == 1 && counter(hh+1)>0
            
            yy = movsum(yEstim,[0 counter(hh+1)])-yEstim;
            yy = yy(1:end-counter(hh+1))/counter(hh+1);
            
        else
            
            yy = yEstim(1+counter(hh+1):end);
            
        end
        
        if isempty(XcondProj)
            
            XX = Xestim(1:end-counter(hh+1),:);
            
        else
            
            XX = [Xestim(1:end-counter(hh+1),:) Xcond(1:end-counter(hh+1),:)];
            
        end
        
        bQRdraws{hhh} = qreg_bayesian( yy, XX,tau,nSave,nBurn,chainStep,bQRoptions);
                bQRtemp = bQRdraws{hhh};
        bQR{hhh} = mean(bQRtemp,3);
        lpQRtemp = zeros(length(tau),nSave/chainStep);
        qXing = zeros(nSave/chainStep,1);
        forProjection = [X(length(y),:) XcondProj];
        
            
        for tt = 1:length(tau)

            lpQRtemp(tt,:) = ...
                forProjection*squeeze(bQRtemp(:,tt,:));

            for dd = 1:nSave/chainStep

                [~,idx] = sort(lpQRtemp(:,dd));

                if ~isequal(idx',1:length(tau))

                    qXing(dd) = 1;

                end

            end

        end
        
        lpQR{hhh} = mean(lpQRtemp,2);   
        lpQRdraws{hhh} = lpQRtemp;
        forProjectionStore{hhh} = forProjection;
        XXstore{hhh} = XX;
        yyStore{hhh} = yy;
        XXlabels{hhh} = [labels labelsTemp];
        
    end
    
end

for hhh = 1:hhL
    
    hh = hhAll(hhh);
    fOutput.(['h' num2str(counter(hh+1))]).bQRdraws = bQRdraws{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).bQR = bQR{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).lpQR = lpQR{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).lpQRdraws = lpQRdraws{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).forProjection = forProjectionStore{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).X = XXstore{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).y = yyStore{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).XXlabels = XXlabels{hhh};
    fOutput.(['h' num2str(counter(hh+1))]).qXing = qXingStore{hhh};
    
    
end


end
