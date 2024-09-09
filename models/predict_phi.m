function phi = predict_phi(phi_d, df, j_d, j, rpm, path_to_model)
% Predict value of phi for given design and operating point
%
% Written by Jordan Eriksen, July 2024
%
% Inputs:
% phi_d         - Design flow coefficient
% df            - Diffusion factor
% j_d           - Design advance ratio
% j             - Operating advance ratio
% rpm           - Operating RPM
% path_to_model - Path to trained model
%
% Outputs:
% phi           - Predicted phi operating point

% Load model into workspace
phi_model = load(path_to_model);

% Prepare table of inputs
T = table(phi_d, df, j_d, j, rpm, 'VariableNames',phi_model.RequiredVariables);

% Predict value
phi = phi_model(T);