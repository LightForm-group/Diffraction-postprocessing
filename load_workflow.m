function [quat_data] = load_workflow(HDF_filepath, phase)

% workflow_name = h5read(HDF_filepath, "/workflow_obj/data/'name'/data");

for num_entries = 1:length(h5info(HDF_filepath, "/element_data/").Groups)
    % list all data as strings
    entry_names{num_entries} = h5info(HDF_filepath, "/element_data/").Groups(num_entries).Name(14:end);

    if contains(entry_names{num_entries}, "volume_element_response") % has VE_response data
        % if entry VE_response get quat data at all incs for task
        quat_data = h5read(HDF_filepath, strcat("/element_data/", entry_names{num_entries}(2:end), ...
            strcat("/data/'phase_data'/data/'",phase,"_orientations'/data/'data'/data/'quaternions'/data")));
    end
end