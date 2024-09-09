function eta_is = predict_eta_is(phi_d, j_d, j, df, phi_op, path_to_model)
% Predict value of eta_is for given design and operating point
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
% eta_is        - Predicted eta_is at operating point

% Load model into workspace
eta_is_model = load(path_to_model);

% Prepare table of inputs
T = table(phi_d, j_d, j, df, phi_op, 'VariableNames',eta_is_model.RequiredVariables);

% Predict value
eta_is = eta_is_model(T);