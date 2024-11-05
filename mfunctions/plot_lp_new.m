function varargout=plot_lp_new(lpQR,y,dates,colour,varargin)

fontSize = 14;

medPos = median(1:size(lpQR,2));
if medPos==floor(medPos)
    fcBands = lpQR(:,[1:medPos-1 medPos+1:end]);
else 
    fcBands = lpQR;
end
fcBands = [repmat(y,1,size(fcBands,2));fcBands];
fcBands(:,2:end) = fcBands(:,2:end)-fcBands(:,1:end-1);

%% Plot area
varargout{1}=figure('units','normalized','outerposition',[0 0.1 .5 0.9]);
h = area(dates',fcBands,min(.995*min(min(fcBands)),1.005*min(min(fcBands))),'LineStyle','none','FaceColor','flat');

% Set colormap
colormap([[1,1,1];fan_chart_color_map(colour,size(fcBands,2)/2)]);

% Plot backdata
hold on
k = plot(dates',[y;NaN(length(dates)-length(y),1)],'Color',colour,'LineWidth',2.5);

axis tight
ylim([-20,20])

yLength = get(gca,'YLim');
line([dates(length(y)),dates(length(y))],yLength,'Color','k')
legend([k h(end-1)],{'Data','Conditional quantiles'},'Location','southoutside','Orientation','horizontal','FontSize',fontSize)

optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont);

if ~isempty(varargin)
   
    SaveFigure(varargin{1},0);

end

end