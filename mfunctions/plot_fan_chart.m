function varargout=plot_fan_chart(fcBands,backData,dates,colour,varargin)

%% Construct date axis
xData = 1:length(dates);

%% Plot area
varargout{1}=figure;
h = area(xData',fcBands,min(.995*min(min(fcBands)),1.005*min(min(fcBands))),'LineStyle','none');

% Set colormap
colormap([[1,1,1];fan_chart_color_map(colour,size(fcBands,2)/2)]);

% Plot backdata
hold on
k = plot(xData,backData,'Color','k','LineWidth',2.5);
% if ~isempty(varargin{1})
%    plot(xData,[prevIR;NaN],'Color','b','LineWidth',2.5);
%    %Will need to take care of additional lines with possibly different lengths
% end

%Change ticks
axis tight

set(gca,'XLim',[min(xData),max(xData)],'YAxisLocation','right')
set(gca,'XTick',min(xData):8:max(xData),'XTickLabel',dates(1:8:end))
legend([k h],'Actual data','Conditional quantiles')

end