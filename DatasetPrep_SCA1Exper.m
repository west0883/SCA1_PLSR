
function [parameters] = DatasetPrep(parameters)

    MessageToUser('Prepping ' , parameters);

    % Pull out data you need
    rest = parameters.rest_data;
    walk = parameters.walk_data;
    
    % if correlations, 
    if strcmp(parameters.data_type, 'correlations')

        % Make dummy variables for rest and walk (response variables). [1 0] =
        % rest. [0 1] = walk.
        responseVariables = [repmat([1 0], size(rest, 3), 1); repmat([0 1], size(walk, 3), 1)];

        % Make explanatory variables (the corrs together)
        explanatoryVariables = cat(3, rest, walk);
    
        % Run Fisher transform
        if isfield(parameters, 'run_FisherTransform') && parameters.run_FisherTransform
    
            explanatoryVariables = atanh(explanatoryVariables);
    
        end 
    
        % reshape explanatoryVarables into vector without repeating values
        explanatoryVariables = [reshape(explanatoryVariables, size(rest, 1) * size(rest,2), size(explanatoryVariables, 3))]';
        indices = find(tril(ones(size(rest, 1)), -1));
        explanatoryVariables = explanatoryVariables(:, indices);
   
    % if fluorescence
    elseif strcmp(parameters.data_type, 'fluorescence')
       
        % Make dummy variables for rest and walk (response variables). [1 0] =
        % rest. [0 1] = walk.
        responseVariables = [repmat([1 0], size(rest, 2), 1); repmat([0 1], size(walk, 2), 1)];
        
        % Make explanatory variables (the fluors together)
        explanatoryVariables = cat(2, rest, walk);

        % transpose explanatory variables
        explanatoryVariables = explanatoryVariables';

    else 
        error('No matching data type.')

    end

    % Remove any rows containing NaNs
    nan_rows = find(any(isnan(explanatoryVariables), 2));
    explanatoryVariables(nan_rows, :) = []; 
    responseVariables(nan_rows, :) = []; 

    % Zscore both variable sets. Keep mu & sigmas for better interprebility
    % of betas later.
    [explanatoryVariables, mu_explanatory, sigma_explanatory] = zscore(explanatoryVariables);
    [responseVariables, mu_response, sigma_response] = zscore(responseVariables);
    dataset.zscoring.explanatoryVariables.mu = mu_explanatory;
    dataset.zscoring.explanatoryVariables.sigma = sigma_explanatory;
    dataset.zscoring.responseVariables.mu = mu_response;
    dataset.zscoring.responseVariables.sigma = sigma_response;

    % For convenience, put both variable sets for this comparison into a stucture for saving. 
    dataset.explanatoryVariables = explanatoryVariables;
    dataset.responseVariables = responseVariables;

    % Put dataset into output structure.
    parameters.dataset = dataset;
    
end