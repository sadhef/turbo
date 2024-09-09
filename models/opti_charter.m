% Assuming 'model_j_adjusted_eta_is' and 'model_phi_gpr' are already loaded and available

%trainedModel = model_j_adjusted_eta_is;
%phi_model = model_phi_gpr;

trainedModel = model_gpr_eta_r2_098;
phi_model = model_phi_gpr;

% Assuming 'model_j_adjusted_eta_is' and 'model_phi_gpr' are already loaded and available

% Parameter ranges
phi_d_range = 0.4:0.05:1;
df_range = 0.25:0.05:0.35;
j_d_range = 0:0.05:1;

% Initialize arrays for storing results
etaIsJ01 = [];
etaIsJ05 = [];
parameterCombinations = [];
allResults = []; % For storing all results

eta_best = 0;
eta_set = [];

tic;
% Iterate through all combinations of parameters for specific j values
for phi_d = phi_d_range
    for df = df_range
        for j_d = j_d_range
            for rpm = [20000]

                if j_d < phi_d * 0.95 

                phi = phi_d ;

                % Evaluate for j=0.1
                j = 0.1;
                T01 = table(phi_d, j_d, df, phi, j, 'VariableNames', trainedModel.RequiredVariables);
                eta_is_01 = trainedModel.predictFcn(T01);
                
                Tp = table(phi_d, df, j_d, j, rpm, 'VariableNames', phi_model.RequiredVariables);
                phi = phi_model.predictFcn(Tp);

                % Evaluate for j=0.5
                j = 0.5;
                T05 = table(phi_d, j_d, df, phi, j, 'VariableNames', trainedModel.RequiredVariables);
                eta_is_05 = trainedModel.predictFcn(T05);
                
                % Store results
                etaIsJ01 = [etaIsJ01; eta_is_01];
                etaIsJ05 = [etaIsJ05; eta_is_05];
                parameterCombinations = [parameterCombinations; [phi_d, df, j_d, 0.1, rpm, phi, eta_is_01]; [phi_d, df, j_d, 0.5, rpm, phi, eta_is_05]];
                allResults = [allResults; phi_d, df, j_d, 0.1, rpm, phi, eta_is_01, 0.5, eta_is_05, eta_is_01 + eta_is_05];

                if (eta_is_01 + eta_is_05) > eta_best
                    eta_best = eta_is_01 + eta_is_05;
                    eta_set = {T01, T05, [eta_is_01, eta_is_05]};
                end
                end
            end
        end
    end
end

% Sort allResults by the combined efficiencies (last column), and select the top 100
[~, sortedIndices] = sort(allResults(:, end), 'descend');
top100Results = allResults(sortedIndices(1:min(100, size(allResults, 1))), :);

% 3D Plot for All Results with phi_d Color Coding
figure;
scatter3(allResults(:, 9), allResults(:, 7), allResults(:, 3), 20, allResults(:, 1), 'filled'); % Change to allResults(:, 1) for phi_d
xlabel('\eta_{is} at j = 0.5');
ylabel('\eta_{is} at j = 0.1');
zlabel('j_d');
title('All Designs with phi_d Color Coding'); % Update title
grid on;
colormap jet; % Keep this for color gradation
colorbar; % Indicates phi_d values the colors represent

% 3D Plot for Top 100 Designs with phi_d Color Coding
figure;
scatter3(top100Results(:, 9), top100Results(:, 7), top100Results(:, 3), 50, top100Results(:, 1), 'filled'); % Change to top100Results(:, 1) for phi_d
xlabel('\eta_{is} at j = 0.5');
ylabel('\eta_{is} at j = 0.1');
zlabel('j_d');
title('Top 100 Designs on the Pareto Front with phi_d Color Coding'); % Update title
grid on;
colormap jet;
colorbar; % Indicates phi_d values the colors represent
xlim([0.6, 0.95]); % Locking X-axis for eta at j = 0.5
ylim([0.6, 0.75]); % Locking Y-axis for eta at j = 0.1

% Assuming the rest of the script is unchanged and top100Results has been correctly calculated

% Extract the relevant information for the top 100 results
% Assuming allResults contains columns in the order: phi_d, df, j_d, j at 0.1, rpm, phi, eta_is_01, j at 0.5, eta_is_05, eta_is_01 + eta_is_05
% For the top 100 results, we want to keep: phi_d, df, j_d, rpm, eta_is_01 at j=0.1, eta_is_05 at j=0.5, combined eta_is
top100ResultsFormatted = top100Results(:, [1, 2, 3, 5, 7, 9, 10]);

% Convert the array to a table for better formatting and saving
colNames = {'Phi_d', 'Df', 'J_d', 'RPM', 'Eta_is_J_0_1', 'Eta_is_J_0_5', 'Combined_Eta_is'};
top100ResultsTable = array2table(top100ResultsFormatted, 'VariableNames', colNames);

% Save the table to an Excel file
filename = 'Top100Designs.xlsx';
writetable(top100ResultsTable, filename);

% Optionally, print the table to MATLAB's command window
disp(top100ResultsTable);

%colNames = {'Phi_d', 'Df', 'J_d', 'RPM', 'Eta_is_J_0_1', 'Eta_is_J_0_5', 'Combined_Eta_is'};
%allResultsNew = array2table(allResults, 'VariableNames', colNames);


filename = 'all_157000_results.xlsx';
writematrix(allResults, filename);

