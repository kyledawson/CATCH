% This script wraps over the GEOS-Chem output and makes the variables in a structure called optics.
% Each file will be saved separately.

filedir = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/datasets/newemissions/';
files = dir([filedir,'*.h5']);

savepath = '../cluster_variables/newemissions/';

for i1 = 1:numel(files)

filename = [filedir,files(i1).name];

optics=get_optics(filename);

sdate = regexp(files(i1).name,'2013\d\d\d\d','match');

savefile = [savepath,'vars_',sdate{:},'.mat'];

save(savefile,'optics','-v7.3')

end

