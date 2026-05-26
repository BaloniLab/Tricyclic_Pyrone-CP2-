% working with iMiceBrain in Matlab 2023a
% load the required toolbox and solver
initCobraToolbox('false');
changeCobraSolver('ibm_cplex', 'all');

%% read iMiceBrain model 
model = readCbModel('/Users/eso1993/Library/CloudStorage/Box-Box/PFOS_Project/drafts_statistics_AND_consensus_draft/models/final_models/iMiceBrain.mat');

% read the file with TPM normalized data 
data = readtable('/Users/eso1993/Library/CloudStorage/Box-Box/CP2_project/data/Stojakovi_paper_young_mice/iMiceBrain_GSE149248_TPM_Normalized.csv', ...
                 'TreatAsEmpty', 'NA', 'ReadVariableNames', true);

for i = 1:width(data)
    if isnumeric(data{:, i})
        nanIdx = isnan(data{:, i});
        data{nanIdx, i} = -1;
    end
end


% set the input parameters 
core = {'BIOMASS_reaction', 'ATPS4mi', 'GLUt6'}; % set the objective function reaction as core

% Extract gene list from the first column
Genes = string(data{:, 1});


%% Loop through each sample (each column except the first)
for i = 2:width(data)
    expressionData.gene = Genes;
    expressionData.value = data{:,i};  % Extract the ith column (expression values)
    
    % Run your mapping function
    [expressionRxns, parsedGPR] = mapExpressionToReactions(model, expressionData);
    
    % fix to set objective function expression to 1 using the rxn index
    expressionRxns(2092,1) = 1;
    expressionRxns(5396,1) = 1;
    expressionRxns(3335,1) = 1;

    
    % Threshold and model construction
    a = median(prctile(expressionRxns, 50));
    tol = 1e-6;
        
    tissueModel = iMAT(model, expressionRxns, 0, a, tol, core);
    tissueModel = updateGenes(tissueModel); % remove inactive genes

    % Save output
    output_folder = '/Users/eso1993/Library/CloudStorage/Box-Box/CP2_project/iMAT_models';
    sampleName = data.Properties.VariableNames{i};
    model_filename = fullfile(output_folder, ['iMAT_' sampleName '.mat']);
    save(model_filename, 'tissueModel');

    fprintf('Processed column: %s\n', sampleName);
end


