function [N, X] = histdiff(Y1, Y2, X)
% HISTDIFF  Efficient histogram of all differences Y1(i) - Y2(j)
%
%   [N, X] = HISTDIFF(Y1, Y2) uses 10 bins.
%   [N, X] = HISTDIFF(Y1, Y2, Nbins) uses Nbins bins.
%   [N, X] = HISTDIFF(Y1, Y2, edges) uses custom edges.
%
%   Note: Unlike MATLAB's HIST, providing a vector X defines bin edges.

    Y1 = Y1(:);
    Y2 = Y2(:);

    % Efficient difference computation: use broadcasting instead of loops
    D = Y1 - Y2';     % Creates matrix of all differences

    if nargin < 3
        % Default: 10 bins
        [N, X] = histcounts(D(:), 10);
        X = (X(1:end-1) + X(2:end))/2;   % Convert edges to centers
    elseif isscalar(X)
        % X = number of bins
        [N, edges] = histcounts(D(:), X);
        X = (edges(1:end-1) + edges(2:end))/2;
    else
        % X = vector of bin edges
        [N, ~] = histcounts(D(:), X);
        X = X(1:end-1);   % Return bin edges as requested by the original function
    end
end
