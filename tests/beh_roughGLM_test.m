% Looping through each of the individual data files
parfor dataFileIdx = 1:length(dataFiles_beh)
    
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_beh{dataFileIdx},dataFileIdx,length(dataFiles_beh))
    
    % We then get the (behavior) filename of the record of interest
    behFilename = [dataFiles_beh{dataFileIdx} '-beh'];
    % and load it into the workspace
    import_data = struct(); import_data = load_behFile(dirs,behFilename);
    
    glm_var = table(); warning off
    
    for trl_i = 10:size(import_data.events.stateFlags_,1)
        glm_var.prev_exp_ssd(trl_i-9) = import_data.events.stateFlags_.LastSsdIdx(trl_i);
        glm_var.prev_rwd(trl_i-9) = import_data.events.stateFlags_.UseRwrdDuration(trl_i-1);
        
        glm_var.curr_rwd(trl_i-9) = import_data.events.stateFlags_.UseRwrdDuration(trl_i);
        
        
        if import_data.events.stateFlags_.IsCancel(trl_i-1) == 1 ||...
                import_data.events.stateFlags_.IsGoCorrect(trl_i-1) == 1
            glm_var.prev_outcome(trl_i-9) = 1;
        elseif import_data.events.stateFlags_.IsCancel(trl_i-1) == 0
            glm_var.prev_outcome(trl_i-9) = 0;
        else
            glm_var.prev_outcome(trl_i-9) = NaN;
        end
        
        if import_data.events.stateFlags_.IsCancel(trl_i-1) == 1 ||...
                import_data.events.stateFlags_.IsGoCorrect(trl_i-1) == 1
            glm_var.curr_outcome(trl_i-9) = 1;
        elseif import_data.events.stateFlags_.IsCancel(trl_i-1) == 0
            glm_var.curr_outcome(trl_i-9) = 0;
        else
            glm_var.curr_outcome(trl_i-9) = NaN;
        end
    end
    
    
    % Fit the GLM
    glm_formula = 'curr_outcome ~ prev_exp_ssd + prev_rwd + curr_rwd + prev_outcome';
    glm_model = [];
    glm_model = fitglm(glm_var, glm_formula, 'Distribution', 'binomial', 'Link', 'logit');
        
    
    glm_out{dataFileIdx} = glm_model;
end




for dataFileIdx = 1:length(dataFiles_beh)
    test(dataFileIdx,:) = glm_out{dataFileIdx}.Coefficients.Estimate';
end
