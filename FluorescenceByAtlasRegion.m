% FLuorescenceByAtlasRegion.m
% Sarah West
%  5/8/23

function [parameters] = FluorescenceByAtlasRegion(parameters)

    MessageToUser('Averaging ', parameters);

    data = parameters.data;
    IC_dim = parameters.IC_dim;
    atlas_ids = parameters.atlas_ids;
    mouse_regions_ordered = parameters.mouse_regions_ordered;

    %%%% average into same atlas region

    % Make a final holding matrix 
    final_timeseries = NaN(size(data, 1), numel(mouse_regions_ordered), size(data, 3));

    % find unique atlas regions
    region_list = unique(atlas_ids);

    % For each atlas region
    for region = region_list

        % get ICs from that region
        these_ICs = find(atlas_ids == region);

        % make an empty list of dimensions 
        C = repmat({':'}, 1, ndims(data));

        % Put IC indices in
        C{IC_dim} = these_ICs;

        % pull out those IC timeseries
        these_timeseries = data(C{:}); 

        % Average across ICs
        this_timeseries = mean(these_timeseries, IC_dim); 

        % Get the final placement of the region
        location = find(mouse_regions_ordered == region);

        % make an empty list of dimensions 
        C = repmat({':'}, 1, ndims(final_timeseries));
        C{IC_dim} = location; 
        final_timeseries(C{:}) = this_timeseries;

    end 
    % Put into output structure
    parameters.data_averaged = final_timeseries;
end 