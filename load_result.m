function [quat_data_all] = load_result(damask_HDF_filepath, phase)

% damask_HDF_filepath = "/Users/user/Downloads/geom_load.hdf5";
damask_result_metadata = h5info(damask_HDF_filepath).Groups;

all_inc_data = [];
% first get list of all increments...
for group_num = 1:length(damask_result_metadata)
    group_name = damask_result_metadata(group_num).Name;
    if contains(group_name, "increment") % has increment data
        inc_num = str2num(group_name(12:end));
        all_inc_data{inc_num+1} = damask_result_metadata(group_num);
    end
end

quat_data_all = [];
% now loop through sorted increments...
for inc_num = 1:length(all_inc_data)
    quat_metadata = all_inc_data{1, inc_num}.Groups(3).Groups.Groups.Datasets(6);
    quat_data_inc = h5read(damask_HDF_filepath, strcat(all_inc_data{inc_num}.Name,"/phase/",phase,"/mechanical/O"));
    quat_data_all = cat(3, quat_data_all, quat_data_inc);
end