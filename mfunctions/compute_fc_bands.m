function varargout=compute_fc_bands(simPaths,simProbsBelowMode)
%Computes the band positions needed to plot a standard fanchart (with 9
%pairs of bands), and optionally also returns the underlying quantiles.
%
%Author: Andrej Sokol

% This function uses a matrix of simulated paths for the variable in
% question, and  vector of probabilities that the variable will be below
% the mode in each quarter, and calculates the matrix that represents the
% positions of each fan chart band.

% A vector from 0.1 to 0.9 in 0.1 steps
pWeights=[1/10:1/10:9/10];

% Calculates the third output of this function. A matrix that contains, for
% each quarter, the probability of the outturn being below each fan band.
% Quarters are rows and fan bands are columns. This is easily possible
% since if, eg, 55% of the probability mass is below the mode and thus 45%
% above, each fan band below the mode will have 5.5% (55%/10) probability
% mass in it and each band above will have 4.5% probability mass in it.
% This is a feature of the two part normal distribution.
varargout{1}=zeros(size(simPaths,2),18);
varargout{1}(:,1:9)=simProbsBelowMode*pWeights;
varargout{1}(:,10:18)=repmat(simProbsBelowMode,1,9)+(1-simProbsBelowMode)*pWeights;
varargout{3} = varargout{1};

% Begins building the second output of this function. A matrix that shows
% the position of each fan chart band for each quarter. The rows are
% quarters and the columns are positions of the bands (or technically, the
% positions of the edges of the bands for each quarter). There are 18
% columns currently. From the output above we know how much probability
% mass (call this 'p') should lie below each fan chart band position. We then use this to
% find the value of the variable such that 'p' fraction of simulations
% ended up being below this value.
for jj=1:size(simPaths,2)
    varargout{1}(jj,:)=quantile(simPaths(:,jj),varargout{1}(jj,:));
end
varargout{2}=varargout{1};

% Calculates the first output of this function: the same as the second but
% in difference form. So if the fan chart band positions for a given
% quarter were [1.1,1.3,1.8,2.0,..] this would transform it into
% [1.1,0.2,0.5,0.2,...] - the first column is the starting the value and
% the rest are the differences
varargout{1}=varargout{1}-[zeros(size(varargout{1},1),1) varargout{1}(:,1:end-1)];
end