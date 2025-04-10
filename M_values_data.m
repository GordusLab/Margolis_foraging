close all;
clear;


% this script runs Model_2025.m twenty times for each value of M, 
% saves and plots the data for the first run and
% calculates the J-S divergence for each. 
%
% JS divergence calculation requires:
%
%data1x = OUTPUT.Exp_Inflection;
%data1y = OUTPUT.Exp_Slopes(1,:) - OUTPUT.Exp_Slopes(2,:);
%data2x = OUTPUT.Model_Inflection;
%data2y = OUTPUT.Model_Slopes(1,:) - OUTPUT.Model_Slopes(2,:);
% 
% Each J-S divergence is stored in a group for each M value 

M = [1000, 100, 10, 1];
num_trials = 20;
js_results = nan(length(M), num_trials);

for i = 1:length(M)
    for j = 1:num_trials
        OUTPUT = Model_2025(M(i))
        
        % Store JS divergence in results matrix

        data1x = OUTPUT.Exp_Inflection;
        data1y = OUTPUT.Exp_Slopes(1,:) - OUTPUT.Exp_Slopes(2,:);
        data2x = OUTPUT.Model_Inflection;
        data2y = OUTPUT.Model_Slopes(1,:) - OUTPUT.Model_Slopes(2,:);
        js_results(i, j) = Jensen_Shannon(data1x,data1y,data2x,data2y);

        % Generate filename e.g., "M1000_Trial1_Model_Output.mat"
        filename = sprintf('M%d_Trial%d_Model_Output.mat', M(i), j);
        filename = fullfile('Model Datasets', filename);
        % Save the OUTPUT structure
        save(filename, 'OUTPUT');
        
        if j == 1 
            Fig2b(OUTPUT, M(i));
            Fig2d(OUTPUT, M(i));
        end
    end
end

save('JS_divergence.mat', 'js_results')
