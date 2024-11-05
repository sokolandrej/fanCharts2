function varargout = plot_qr_coefficients(bQR,lb,ub,bOLS,lbOLS,ubOLS,quantiles,ci,labels,colour,varargin)
% Plots regression coefficients by quantile.

sb = size(bQR,1);
pp = numSubplots(size(bQR,1));
fontSize = 14;

if isempty(labels)
    
   labels = cell(sb,1); 
   
   for ii = 1:sb 
      
       labels{ii} = ['b' num2str(ii)];
       
   end
    
end

if sb == 1
    
    varargout{1} = figure('units','normalized','outerposition',[0.1 0.1 .45 0.9]);

else
    
    varargout{1} = figure('units','normalized','outerposition',[0.1 0.1 .9 0.9]);

end

for ii = 1:sb
    
    subplot(pp(1),pp(2),ii)
   
    H = PlotSwathe(bQR(ii,:)', [ub(ii,:)' lb(ii,:)'], rgb(colour));hold on
    handle(1)=H.bar; handle(2) = H.patch; 
% if size(bOLS,2) == 1
%     handle(3) = plot(bOLS(ii)*ones(length(quantiles),1),'LineWidth',0.5,'LineStyle','-','Color',rgb('black')); hold on
%     plot(lbOLS(ii)*ones(length(quantiles),1),'LineWidth',0.5,'LineStyle','--','Color',rgb('black')); hold on
%     plot(ubOLS(ii)*ones(length(quantiles),1),'LineWidth',0.5,'LineStyle','--','Color',rgb('black')); hold on
% else
%     handle(3) = plot(bOLS(ii,:),'LineWidth',2,'LineStyle','--','Color',rgb('black')); hold on
%     plot(lbOLS(ii,:),'LineWidth',2,'LineStyle','--','Color',rgb('black')); hold on
%     plot(ubOLS(ii,:),'LineWidth',2,'LineStyle','--','Color',rgb('black')); hold on
% end

plot(zeros(length(quantiles),1),'LineWidth',0.5,'LineStyle','-','Color',rgb('red'));
    title(labels(ii)) 

%     set(gca,'xTick',[1 fix(median(1:length(quantiles))) length(quantiles)],...
%         'xTickLabel',quantiles([1 fix(median(1:length(quantiles))) length(quantiles)]),...
%         'xLim',[0 length(quantiles)+1],'Layer','top');
set(gca,'xTick',[1 round([length(quantiles)/4 length(quantiles)/2 3*length(quantiles)/4 length(quantiles)])],...
        'xTickLabel',quantiles([1 round([length(quantiles)/4 length(quantiles)/2 3*length(quantiles)/4 length(quantiles)])]),...
        'xLim',[0 length(quantiles)+1],'yTick',[-2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5],'yTickLabel',[-2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5],'yLim',[-2.5 1.5],'Layer','top');

axis tight  
% ylim([-2.5 1.5])

set(gca, 'XLimSpec', 'Tight');
    optfont = FigFontOption(fontSize) ; FigFont(optfont); %optfont.title_weight = 'bold'
         xlabel('\tau')
%     ylabel('Percentage points')
end
% opt=LegOption; opt.handle = handle; LegSubplot({'Point estimate',[num2str(ci) '% C.I.'],'OLS estimate'},opt);
% varargout{2} = opt;

   
if ~isempty(varargin)
   
    SaveFigure(varargin{1},0);

end

end


