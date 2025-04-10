function js_divergence = Jensen_Shannon(data1x,data1y,data2x,data2y)


% data1 is reference

x_bins = freedman_diaconis(data1x);
y_bins = freedman_diaconis(data1y);

P = histcounts2(data1x, data1y, x_bins,y_bins,'Normalization', 'probability');
Q = histcounts2(data2x, data2y, x_bins,y_bins,'Normalization', 'probability');

% Calculate J-S

M = 0.5 * (P + Q);

% Calculate KL divergence D_KL(P || M)
D_KL_P_M = sum(P .* log(P ./ M), 'all', 'omitnan');  % 'omitnan' to handle zero values

% Calculate KL divergence D_KL(Q || M)
D_KL_Q_M = sum(Q .* log(Q ./ M), 'all', 'omitnan');  % 'omitnan' to handle zero values

% Jensen-Shannon Divergence: 0.5 * (D_KL(P || M) + D_KL(Q || M))
js_divergence = 0.5 * (D_KL_P_M + D_KL_Q_M);


end

function binEdges = freedman_diaconis(data)


Q1 = prctile(data, 25);  % 25th percentile
Q3 = prctile(data, 75);  % 75th percentile
IQR = Q3 - Q1;           % Interquartile range

n = length(data);          % Number of data points
binWidth = (2 * IQR) / (n^(1/3));  % Freedman-Diaconis bin width

binEdges = min(data):binWidth:max(data);  % Generate bin edges

binEdges = [-Inf, binEdges, Inf];

end