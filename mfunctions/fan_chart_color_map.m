function colorMap=fan_chart_color_map(baseColor,nBands)
%This function constructs a color map for fan charts and zigurrat charts
%from an RGB input color (note that RGB values must lie between 0 and 1,
%which means that you might divide your triplet by 255).
%
%Author: Andrej Sokol

% if mod(nBands,2)==0
%     error('The number of bands must be odd');
%     return
% end
idx=nBands-1:-1:1;
col1=linspace(.925,baseColor(1),nBands);
col1=[col1 col1(idx)];
col2=linspace(.925,baseColor(2),nBands);
col2=[col2 col2(idx)];
col3=linspace(.925,baseColor(3),nBands);
col3=[col3 col3(idx)];
colorMap=[col1' col2' col3'];

end