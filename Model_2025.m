function OUTPUT = Model_2025(M0)

% This function models reorientation data based on a simple two variable
% ODE model, where Reorientation Rate (dRdt) is a function of a decaying
% reorientation drive signal (Y) in the following manner:
%
% dRdt = alpha*Y
% dYdt = -gamma*Y
% The boundary condition for Y is that Y(t=0) is drawn from a lognormal
% distribution with mean = 1 & std = 0.5. This is done so that the model
% variance matches the experimental variance. Parameters alpha and
% gamma are optimized so that the residual of the average reorientation
% data and average model data is minimized.
% 
% Data are subsequently compared using findchangepts, where the difference
% in slopes before and after the inflection point are compared.
%
% INPUT
% data: cumulative sum of reorientation
%       rows: frames (3fps, 8100 frames= 45 min)
%       columns: individuals
%
% OUTPUT
%
% OUTPUT.Exp_Data = Experimental Data;
% OUTPUT.Time = Experimental Time (3 fps);
% OUTPUT.Model_Parameters = Optimized alpha, beta & gamma.
% OUTPUT.PredictedRate = Predicted reorientation rate based on optimized
% parameters.
% OUTPUT.Model_Data = Discrete reorientation data modelled with the
% Gillespie algorithm, using the optimized parameters.
% OUTPUT.Exp_Inflection = Experimental Inflection points using
% findchangepts.
% OUTPUT.Exp_Slopes = Experimental Slopes using findchangepts. The first
% row is the slope before the inflection, and the second row is the slope
% after the inflection.
% OUTPUT.Model_Inflection = Model Inflection points using findchangepts.
% OUTPUT.Model_Slopes = Model Slopes using findchangepts. The first row is 
% the slope before the inflection, and the second row is the slope after 
% the inflection.
% OUTPUT.MODEL_XCORR_INDICES = The indices of model data that most highly
% correlate with experimental data. Each index of this vector represent
% each index of the experimental data. Each value within the vector is the
% index for the model data that most closely correlates with that vector
% index.
%
% Andrew Margolis & Andrew Gordus
%
% August 2023

load('N2_cumevents_matrix.mat');

data = N2_cumevent_matrix;
clear N2_cumevent_matrix;

% Clear out first 2 minutes as recommended by Lopez-Cruz

data(1:360,:) = [];
[d1,d2] = size(data);
data = data - ones(d1,1)*data(1,:); % resets count to start after 2 min.
time = (0:size(data,1)-1)/(3*60); % time in minutes

OUTPUT.Exp_Data = data;
OUTPUT.Time = time;

%----------------------%
% PARAMETER ESTIMATION %
%----------------------%

% Produce average reorientation rate using a moving 2 minute window

[mean_rate, initial_rates] = moving_average(OUTPUT);
% 

% Initial test rates. All rates are min^-1;
alpha = 1.5;
gamma = 0.1;
beta = mean_rate(end);
parameters = [alpha,gamma,beta];

% M0 is M at time = zero.


% Fit model parameters to mean rate of experimental data

t = transpose(OUTPUT.Time(361:end) - OUTPUT.Time(361));
y = transpose(mean_rate);

fitparameters = nlinfit(t,y,@ratemodel,parameters);

OUTPUT.PredictedRate = ratemodel(fitparameters,time);

OUTPUT.Model_Parameters = fitparameters;

% Generate d2 number of in silico reorientation rates
% Initial rates drawn from distribution of initial rates from experimental
% data.

[Time,Reorientations,M] = gillespie_search(fitparameters,time(end),d2,initial_rates,M0);

% Grab inter-orientation time intervals for model

OUTPUT.model_time_intervals = [];
for j=1:length(Reorientations)
    reos = Reorientations{j};
    reos_diff = [0;diff(reos)];
    temp_t = Time{j};
    reorient_times = temp_t(reos_diff==1);
    OUTPUT.model_time_intervals = [OUTPUT.model_time_intervals; diff(reorient_times)];
end

OUTPUT.Model_Data_Raw = Reorientations;
OUTPUT.Model_Data_Raw_Time = Time;
OUTPUT.M_raw = M;

% Grab inter-orientation time intervals for experimental data

OUTPUT.exp_time_intervals = [];

for j=1:size(OUTPUT.Exp_Data,2)
    reos = OUTPUT.Exp_Data(:,j);
    reos_diff = [0;diff(reos)];
    reorient_times = OUTPUT.Time(reos_diff==1);
    OUTPUT.exp_time_intervals = [OUTPUT.exp_time_intervals, diff(reorient_times)];
end
    
% Rescale model time to experimental time
Reorientations = fixtime(Time,Reorientations,time);
M = fixtime(Time,M,time);

OUTPUT.Model_Data = Reorientations;
OUTPUT.M_data = M;

clear Time Reorientations M


%-----------------------------------------------%
% FIND INFLECTION POINTS & CORRESPONDING SLOPES %
%-----------------------------------------------%

% findchangepts 

OUTPUT.Exp_Inflection = NaN(1,d2);
OUTPUT.Exp_Slopes = NaN(2,d2);

OUTPUT.Model_Inflection = NaN(1,d2);
OUTPUT.Model_Slopes = NaN(2,d2);

for m=1:d2
    
    if sum(isnan(data(:,m))) < d1
        ipt = findchangepts(data(:,m));
        OUTPUT.Exp_Inflection(m) = ipt;
        
        b1 = [ones(ipt,1),(0:ipt-1)']\data(1:ipt,m);
        b2 = [ones(d1-ipt+1,1),(0:d1-ipt)']\(data(ipt:d1,m) - data(ipt,m));
        OUTPUT.Exp_Slopes(1,m) = 3*60*b1(2); % slopes min^-1
        OUTPUT.Exp_Slopes(2,m) = 3*60*b2(2);
    else
        OUTPUT.Exp_Inflection(m) = 1;
        OUTPUT.Exp_Slopes(1,m) = 0;
        OUTPUT.Exp_Slopes(2,m) = 0;
    end
    
    if sum(isnan(OUTPUT.Model_Data(:,m))) < d1
        ipt = findchangepts(OUTPUT.Model_Data(:,m));
        OUTPUT.Model_Inflection(m) = ipt;
        
        b1 = [ones(ipt,1),(0:ipt-1)']\OUTPUT.Model_Data(1:ipt,m);
        b2 = [ones(d1-ipt+1,1),(0:d1-ipt)']\(OUTPUT.Model_Data(ipt:d1,m) - OUTPUT.Model_Data(ipt,m));
        OUTPUT.Model_Slopes(1,m) = 3*60*b1(2); % slopes min^-1
        OUTPUT.Model_Slopes(2,m) = 3*60*b2(2);
    else
        OUTPUT.Model_Inflection(m) = 1;
        OUTPUT.Model_Slopes(1,m) = 0;
        OUTPUT.Model_Slopes(2,m) = 0;
    end
end


end



%------------------------%
% MOVING AVERAGE OF RATE %
%------------------------%

function [mean_rate, initial_rates] = moving_average(OUTPUT)

window = 360;
length = size(OUTPUT.Time,2) - window;

mean_rate = zeros(1,length);

count = 1;
for j=361:size(OUTPUT.Time,2)

    start = j - window;
    mean_rate(count) = mean((OUTPUT.Exp_Data(j,:) - OUTPUT.Exp_Data(start,:)) ./2 );
    count = count+1;

    if j==361
        initial_rates = (OUTPUT.Exp_Data(j,:) - OUTPUT.Exp_Data(start,:)) ./2 ;
    end
end

end

%-------%
% MODEL %
%-------%
function Reorientation_rate = ratemodel(parameters,time)

alpha = parameters(1);
gamma = parameters(2);
beta = parameters(3);

Reorientation_rate = alpha*exp(-gamma*time) + beta;

end

%-----------%
% GILLESPIE %
%-----------%

function [Time,Reorientations,M] = gillespie_search(parameters,time_limit,number,initial_rates,M0)

% Reactions
%
% R -> R+1, rate = (alpha*y), R=R+1
% Y -> null, rate = -gamma*Y, Y = Y-1, this is loss of signal that promotes
% local search
% X = [R, Y]
% Xstart = [#Reorientations #Modulators];
%
% I chose #Modulators to be a large aribtrary number so that the decay is fairly
% continuous.

% number = length(initial_rates);
Xstart = [0 M0];

% Initial alpha value weighted %

m = length(initial_rates);
random_index = randi(m,1,m);
Yr = initial_rates(random_index)./mean(initial_rates);

if M0==1
    Yr(Yr > 1) = 1;
    Yr(Yr < 1) = 0;
end

Time = cell(number,1);
Reorientations = cell(number,1);
M = cell(number,1);


stoich_matrix = [1 0;0 -1]; % This is the stochiometry of reaction R & Y

alpha = parameters(1)/Xstart(2);
% gamma = parameters(2);
gamma = 0.5;
beta = parameters(3);

maxlength = 10000;


for cycle =1:number
    T = NaN(maxlength,1);
    X = NaN(maxlength,length(Xstart));
    
    T(1) = 0;
    
    if M0==1
        X(1,:) = Xstart;
    else
        X(1,:) = [Xstart(1) ceil(Xstart(2)*Yr(cycle))];
    end
    
    count = 1;
    while T(count,1) < time_limit % Do the simulation for 50 min window

        % probability of different reactions

        a = [(alpha*X(count,2)) + beta,gamma*X(count,2)];
        a0 = sum(a);
        r = rand(1,2);
        
        % step length
        tau = -log(r(1))/a0;
        
        %reaction choice
        mu = find((cumsum(a) >= r(2)*a0), 1,'first');
        
        % Update time and reactants
        T(count+1) = T(count)+tau;
        X(count+1,:) = X(count,:) + stoich_matrix(mu,:);
        count = count+1;
        if count==maxlength
            disp(['Reached maximum length of ',num2str(maxlength)]);
            break
        end
    end
    
    T(count:end) = [];
    X(count:end,:) = [];
    Time(cycle) = {T};
    Reorientations(cycle) = {X(:,1)};
    M(cycle) = {X(:,2)};
end

end

%---------------------------------------%
% PUT MODEL TIME INTO EXPERIMENTAL TIME %
%---------------------------------------%

function NewReorientations = fixtime(Model_Time,Reorientations,Exp_Time)

NewReorientations = zeros(length(Exp_Time),length(Reorientations));

for j=1:length(Model_Time)
    reos = Reorientations{j};
    for m=2:length(Exp_Time)
        n = (Exp_Time(m-1)<Model_Time{j}) & (Model_Time{j}<Exp_Time(m));
        if ~isempty(max(reos(n)))
            NewReorientations(m,j) = max(reos(n));
        else
            NewReorientations(m,j) = NewReorientations(m-1,j);
        end
    end
    disp(['Rescaling ',num2str(j),' of ',num2str(length(Model_Time))]);
end

end
