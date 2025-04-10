function output = Fig2d(OUTPUT, M)


% This code takes experimental and model data from the output of Model.m 
% and plots individual worm traces, the distribution of traces based on 
% inflection point and change in slope, the aggregate trace, and 
% histograms for the distribution of inflection points and slope
% differences.
% 
% Individual worm traces are plotted from each quadrant of the 
% two-dimensional distribution.


%% Figures 1a,1b,1c, and (All individual experimental traces)
close all;

% Organize traces based on inflection point and pull out traces around
% quadrant centroids (Xcentroid was set to different values to pull 
% from around other points in the distribution)

trace_coord = cat(2,OUTPUT.Time(OUTPUT.Exp_Inflection).',OUTPUT.Exp_Slopes(1,:).' - OUTPUT.Exp_Slopes(2,:).');
slopes_median = median(OUTPUT.Exp_Slopes(1,:) - OUTPUT.Exp_Slopes(2,:));
inflection_median = median(OUTPUT.Time(OUTPUT.Exp_Inflection));

% Divide traces into quadrants and find traces 
[row, col] = find(trace_coord(:,2)> slopes_median);
AB_coord = trace_coord(row,:);
[row, col] = find(trace_coord(:,2)< slopes_median);
CD_coord = trace_coord(row,:);
[rowA, colA] = find(AB_coord(:,1) < inflection_median); 
Acoord = AB_coord(rowA,:); 
Acentroid = [mean(Acoord(:,1))  mean(Acoord(:,2))];
A_Exp_Traceidx = dsearchn(trace_coord,Acentroid); % determine centroid
A_Exp_TraceCoord = trace_coord(A_Exp_Traceidx,:);
[rowB, colB] = find(AB_coord(:,1) > inflection_median);
Bcoord = AB_coord(rowB,:);
Bcentroid = [mean(Bcoord(:,1))  mean(Bcoord(:,2))];
B_Exp_Traceidx = dsearchn(trace_coord,Bcentroid);
B_Exp_TraceCoord = trace_coord(B_Exp_Traceidx,:); % determine closest datapoint to centroid
[rowC, colC] = find(CD_coord(:,1) < inflection_median);
Ccoord = CD_coord(rowC,:);
Ccentroid = [mean(Ccoord(:,1))  mean(Ccoord(:,2))];
C_Exp_Traceidx = dsearchn(trace_coord,Ccentroid);
C_Exp_TraceCoord = trace_coord(C_Exp_Traceidx,:);
[rowD, colD] = find(CD_coord(:,1) > inflection_median);
Dcoord = CD_coord(rowD,:);
Dcentroid = [mean(Dcoord(:,1))  mean(Dcoord(:,2))];
D_Exp_Traceidx = dsearchn(trace_coord,Dcentroid);%19, 0
D_Exp_TraceCoord = trace_coord(D_Exp_Traceidx,:);

% Plot traces
figure(1);
subplot(2,2,1)
trace = [OUTPUT.Time.' , (OUTPUT.Exp_Data(:,A_Exp_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['A' '(' join(string(A_Exp_TraceCoord),' , '), ')'],' '));

subplot(2,2,2)
trace = [OUTPUT.Time.' , (OUTPUT.Exp_Data(:,B_Exp_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['B' '(' join(string(B_Exp_TraceCoord),' , '), ')'],' '));            

subplot(2,2,3)
trace = [OUTPUT.Time.' , (OUTPUT.Exp_Data(:,C_Exp_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['C' '(' join(string(C_Exp_TraceCoord),' , '), ')'],' '));            

subplot(2,2,4)
trace = [OUTPUT.Time.' , (OUTPUT.Exp_Data(:,D_Exp_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['D' '(' join(string(D_Exp_TraceCoord),' , '), ')'],' '));            

filename = sprintf('M%d_Fig2d_ExpTraces.fig', M);
foldername = sprintf('M%d', M);
filename = fullfile(strcat('Figures/',foldername), filename);
saveas(gcf, filename);
close(gcf);


%% Figures 1d and 1f (All individual model traces inset with distribution)

trace_coord = cat(2,OUTPUT.Time(OUTPUT.Model_Inflection).',OUTPUT.Model_Slopes(1,:).' - OUTPUT.Exp_Slopes(2,:).');

% Divide traces into quadrants and find traces around coordinates of
% experimental traces

A_Model_Traceidx = dsearchn(trace_coord,A_Exp_TraceCoord);
A_Model_TraceCoord = trace_coord(A_Model_Traceidx,:);
B_Model_Traceidx = dsearchn(trace_coord,B_Exp_TraceCoord);
B_Model_TraceCoord = trace_coord(B_Model_Traceidx,:);
C_Model_Traceidx = dsearchn(trace_coord,C_Exp_TraceCoord);
C_Model_TraceCoord = trace_coord(C_Model_Traceidx,:);
D_Model_Traceidx = dsearchn(trace_coord,D_Exp_TraceCoord);
D_Model_TraceCoord = trace_coord(D_Model_Traceidx,:);

% Plot traces
figure(2);
subplot(2,2,1)

trace = [OUTPUT.Time.' , (OUTPUT.Model_Data(:,A_Model_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['A' '(' join(string(A_Model_TraceCoord),' , '), ')'],' '));

subplot(2,2,2)
trace = [OUTPUT.Time.' , (OUTPUT.Model_Data(:,B_Model_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['B' '(' join(string(B_Model_TraceCoord),' , '), ')'],' '));

subplot(2,2,3)
trace = [OUTPUT.Time.' , (OUTPUT.Model_Data(:,C_Model_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['C' '(' join(string(C_Model_TraceCoord),' , '), ')'],' '));

subplot(2,2,4)
trace = [OUTPUT.Time.' , (OUTPUT.Model_Data(:,D_Model_Traceidx))];
newtrace = figure_trace(trace);
scatter(newtrace(:,1),newtrace(:,2),20);
title(join(['D' '(' join(string(D_Model_TraceCoord),' , '), ')'],' '));

sgtitle('Model')

filename = sprintf('M%d_Fig2d_ModelTraces.fig', M);
foldername = sprintf('M%d', M);
filename = fullfile(strcat('Figures/',foldername), filename);
saveas(gcf, filename);
close(gcf);


% fig = gcf;
% fig.PaperPositionMode = 'auto';
% saveas(gcf,'Model Traces','svg')

%% Figures 1e and 1g (Cumulative Trace)

figure;

% Cumulative Trace
subplot(2,2,1)
h1 = plot(OUTPUT.Time,nanmean(OUTPUT.Exp_Data,2),'b'); hold on;
plot(OUTPUT.Time,nanmean(OUTPUT.Exp_Data,2)+nanstd(OUTPUT.Exp_Data,0,2),'b--');
plot(OUTPUT.Time,nanmean(OUTPUT.Exp_Data,2)-nanstd(OUTPUT.Exp_Data,0,2),'b--');
h2 = plot(OUTPUT.Time,nanmean(OUTPUT.Model_Data,2),'r');
plot(OUTPUT.Time,nanmean(OUTPUT.Model_Data,2)+nanstd(OUTPUT.Model_Data,0,2),'r--');
plot(OUTPUT.Time,nanmean(OUTPUT.Model_Data,2)-nanstd(OUTPUT.Model_Data,0,2),'r--');
ylim([0,40]);
xlabel('Time (minutes)');
ylabel('Number of Reorientations');
title('Cumulative Reorientations');
legend([h1,h2],{'Exp Data','Model Data'});

%% Figures 1d and 1f (Scatterplot with Contours)

% Scatterplot
subplot(2,2,2)
h1 = plot(OUTPUT.Time(OUTPUT.Exp_Inflection),OUTPUT.Exp_Slopes(1,:) - OUTPUT.Exp_Slopes(2,:),'b.','MarkerSize',0.5,'DisplayName','Experimental Data');
hold on;
h2 = plot(OUTPUT.Time(OUTPUT.Model_Inflection),OUTPUT.Model_Slopes(1,:) - OUTPUT.Model_Slopes(2,:),'r.','MarkerSize',0.5,'DisplayName','Model Data');
xl = xlim;
yl = ylim;

% Denote coordinates of sample traces
plot(A_Model_TraceCoord(1),A_Model_TraceCoord(2),'k.');
plot(A_Exp_TraceCoord(1),A_Exp_TraceCoord(2),'k.');
plot(D_Model_TraceCoord(1),D_Model_TraceCoord(2),'k.');
plot(D_Exp_TraceCoord(1),D_Exp_TraceCoord(2),'k.');

% % Experimental contour
% [N,C] = hist3([OUTPUT.Time(OUTPUT.Exp_Inflection).',(OUTPUT.Exp_Slopes(1,:) - OUTPUT.Exp_Slopes(2,:)).'],[15,15]);
% Nnormalized = N / sum(N,'all');
% [Mexp,cexp] = contour(C{1},C{2},Nnormalized.',3,'LineColor','b','LineWidth',0.75,'ShowText','on');
% 
% % Model contour
% [N,C] = hist3([OUTPUT.Time(OUTPUT.Model_Inflection).',(OUTPUT.Model_Slopes(1,:) - OUTPUT.Model_Slopes(2,:)).'],[15,15]);
% Nnormalized = N / sum(N,'all');
% [Mmod,cmod] = contour(C{1},C{2},Nnormalized.',3,'LineColor','r','LineWidth',0.75,'ShowText','on');

xlabel('Inflection Point (minutes)');
ylabel('Slope 1 - Slope 2 (min^-^1)');
title('Relative Change in Slope at Inflection Points');
legend([h1,h2],{'Exp Data','Model Data'});

%% Figures 1d and 1f (Histograms)

% Change in Slope Histogram - Experimental
subplot(2,2,3)
edges = linspace(yl(1),yl(2));
h = histogram(OUTPUT.Exp_Slopes(1,:) - OUTPUT.Exp_Slopes(2,:),edges,'FaceColor','b','Normalization','probability');
hold on

% Fit to log-normal curve
Exp_PositiveDeltaSlopes = OUTPUT.Exp_Slopes(1,:) - OUTPUT.Exp_Slopes(2,:);
Exp_PositiveDeltaSlopes(find(Exp_PositiveDeltaSlopes<=0))=[];
OUTPUT.ExpPar_logn = lognfit(Exp_PositiveDeltaSlopes);
area = sum(h.Values)*h.BinWidth;
% plot(edges,area*lognpdf(edges,OUTPUT.ExpPar_logn(1),OUTPUT.ExpPar_logn(2)),'b','LineWidth',0.5)
hold on
xlabel('Slope 1 - Slope 2 (min^-^1)');
ylabel('Frequency');

% Change in Slope Histogram - Model
h = histogram(OUTPUT.Model_Slopes(1,:) - OUTPUT.Model_Slopes(2,:),edges,'FaceColor','r','Normalization','probability');
hold on

% Fit to log-normal curve
Model_PositiveDeltaSlopes = OUTPUT.Model_Slopes(1,:) - OUTPUT.Model_Slopes(2,:);
Model_PositiveDeltaSlopes(find(Model_PositiveDeltaSlopes<=0))=[];
OUTPUT.ModelPar_logn = lognfit(Model_PositiveDeltaSlopes);
area = sum(h.Values)*h.BinWidth;
% plot(edges,area*lognpdf(edges,OUTPUT.ModelPar_logn(1),OUTPUT.ModelPar_logn(2)),'r','LineWidth',0.5)

% Inflection Points Histogram - Experimental
subplot(2,2,4)
edges = linspace(xl(1),xl(2));
h = histogram(OUTPUT.Time(OUTPUT.Exp_Inflection),edges,'FaceColor','b','Normalization','probability');
hold on

% Fit to Gaussian curve
% f = fit(h.BinEdges(2:end).',h.Values.','gauss1');
% plot(f,'b');
hold on

% Inflection Points Histogram - Model
h = histogram(OUTPUT.Time(OUTPUT.Model_Inflection),edges,'FaceColor','r','Normalization','probability');

% Fit to Gaussian curve
% f = fit(h.BinEdges(2:end).',h.Values.','gauss1');
% plot(f,'r');

xlabel('Inflection Point (minutes)');
ylabel('Frequency');

filename = sprintf('M%d_Fig2d_Scatter.fig', M);
foldername = sprintf('M%d', M);
filename = fullfile(strcat('Figures/',foldername), filename);
saveas(gcf, filename);
close(gcf);


% fig = gcf;
% fig.PaperPositionMode = 'auto';
% saveas(gcf,'Scatterplot_CumulativeTrace_Histograms','svg')


% This function removes redundant points from traces so that the traces
% only contain timepoints where the number of cumulative reorientations
% changes

function newTrace = figure_trace(trace)

% Initialize variables
newTrace = trace(1,:);
prev_y = trace(1,2);

% Loop through each row of the trace array
for i = 2:size(trace,1)
    % Check if the y-value changes from the previous value
    if trace(i,2) ~= prev_y
        % If it does, add the current row to the newTrace array
        newTrace = [newTrace; trace(i,:)];
        prev_y = trace(i,2);
    end
end


