function psi = predict_psi(phi_d, j_d, j, df, phi_op, path_to_model)
% Predict value of psi for given design and operating point
%
% Written by Jordan Eriksen, July 2024
%
% Inputs:
% phi_d         - Design flow coefficient
% df            - Diffusion factor
% j_d           - Design advance ratio
% j             - Operating advance ratio
% phi_op        - Operating flow coefficient
% path_to_model - Path to trained model
%
% Outputs:
% psi           - Predicted psi operating point

% Load model into workspace
psi_model = load(path_to_model);

% Prepare table of inputs
T = table(phi_d, j_d, j, df, phi_op, 'VariableNames',psi_model.RequiredVariables);

% Predict value
psi = psi_model(T);