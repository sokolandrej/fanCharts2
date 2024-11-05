function [y,x,x_f] = realtime_XY(data,vintage,nfore,lags)


% Transform to stationarity
all_series = [];
for i = 1:size(data,1)
    tcode = data{i,3};
    if tcode == 5    % Annualized growth rates
        all_series(:,i) = 400*transx(data{i,1}(1:vintage+111,vintage),tcode);
    elseif tcode == 4 || tcode == 2 % Logs or differences
        all_series(:,i) = 4*transx(data{i,1}(1:vintage+111,vintage),tcode);
    else             % Levels
        all_series(:,i) = data{i,1}(1:vintage+111,vintage);
    end
end

% Specify which series are LHS and which are RHS
Yraw_full = yfcsta(all_series(:,1),data{1,3}(1),nfore);
Yraw_full_lags = all_series(:,1);
Zraw = all_series(:,2:end);

% Now correct observations after taking lagged values
Yraw_full = Yraw_full(1:end-nfore,:);
% Keep for forecasting
Yraw_full_lags_fore = [Yraw_full_lags mlag2(Yraw_full_lags,lags-1)];
Yraw_full_lags_fore = Yraw_full_lags_fore(end,:);
Zfore = Zraw(end,:);
% Keep for estimation
Yraw_full_lags = Yraw_full_lags(1:end-nfore,:);
Zraw = Zraw(1:end-nfore,:);

% Create RHS variables, x
ylag = mlag2(Yraw_full_lags,lags-1);  % This routine creates the lagged dep variables   
if lags>0
    Xraw = [ones(size(ylag(lags:end,:),1),1) Yraw_full_lags(lags:end,:) ylag(lags:end,:) Zraw(lags:end,:)];
    Xfore = [1, Yraw_full_lags_fore, Zfore];
    Yraw = Yraw_full(lags:end,:);
elseif lags==0
    Xraw = [ones(size(ylag(1:end,:),1),1) Zraw(1:end,:)];
    Xfore = [1, Zfore];
    Yraw = Yraw_full(1:end,:);
end

y = Yraw;
x = Xraw;
x_f = Xfore;