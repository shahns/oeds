function plot_cov_ellipse(mu, covMat, varargin)
% PLOT_COV_ELLIPSE plots covariance matrix as a linear a transformation
%        specified by mean M and covariance C Dimentionality of
%        should match M and C
%
% Syntax:   PLOT_COV_ELLIPSE(mu, covMat, 'showMean' true, ...
%               'legendText', "abc", 'labels', ["x","y"])
%
% Inputs:
%   REQUIRED:
%    mu - mean of the distribution
%    covMat - covariance of the distribution
%
%   OPTIONAL PARAMETERS:
%    'showMean' - BOOLEAN to show mean
%    'legendText' - STRING ARRAY to show legend
%    'labels' - STRING ARRAY to show axis labels x followed by y
%               default valye ["x", "y"]
%    'color' - specify color for ellipse only std matlab colors
%                     are valid
defaultShowMean = false;
defaultLegendText = '';
defaultLabels = ["x", "y"];
defaultColor = 'm';
validateLabels = @(x) isequal(size(x), [1, 2]);
validateColor = @(x) isequal(size(x), [1, 1])&& ischar(x)|| isequal(size(x), [1, 3]);
p = inputParser;
addRequired(p, 'mu');
addRequired(p, 'covMat');
addParameter(p, 'showMean', defaultShowMean, @islogical);
addParameter(p, 'legendText', defaultLegendText);
addParameter(p, 'labels', defaultLabels, validateLabels);
addParameter(p, 'color', defaultColor, validateColor);
parse(p, mu, covMat, varargin{:});
if ~(all(numel(p.Results.mu) == size(p.Results.covMat)))
    error('Dimensionality of mu and covMat must match');
end
thetas = linspace(0, 360, 720);
[V, D] = eig(p.Results.covMat);
a = sqrt(diag(D));
points = V * [cosd(thetas) * a(1); sind(thetas) * a(2)] + p.Results.mu;
plot(points(1, :), points(2, :), p.Results.color, 'LineWidth', 1);
hold on
if (p.Results.showMean)
    plot(mu(1), mu(2), '*', 'MarkerSize', 5, 'color' , p.Results.color);
end
if ~(strcmp(p.Results.legendText, ''))
    legend(p.Results.legendText);
end
xlabel(p.Results.labels(1));
ylabel(p.Results.labels(2));
end
