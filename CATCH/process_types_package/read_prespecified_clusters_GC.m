function [means,cov,means_column_info] = read_prespecified_clusters_GC(version)
% This function will read in Sharon's cluster data for the centroids and
% ellipses of the mahalanobis distance metric.
% set the version number, 1 = oldest
%                         2 = one with distribution data
%                         3 = selected file rawdata and log of it.

if ~isstr(version)||isnumeric(version)
 version = num2str(version);
end

path = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/process_types_package/';
path2 = ['Aerosol_Clusters_v',...
    sprintf('%s',version),'.txt'];
[~,svers] = fileparts(path2);
fprintf('Using version: %s\n',svers);
load([path,'colormaps.mat'])
Class{3}='DustyMix';
Class{7}='FreshSmoke';
Class{8}='PollutedMaritime';
Class{9}='PureDust';

nheaders = 4;
nvars = 5;
delim = ',';
frmt = repmat('%f',1,nvars);

fid = fopen([path,path2],'r');
means = textscan(fid,['%d',frmt,'\n'],5,'headerlines',nheaders,'delimiter',delim,...
    'collectoutput',1);
types = Class(means{1}); means = means{2};

for i1 = 1:numel(types)
    temp = textscan(fid,frmt,nvars,'headerlines',1,'delimiter',delim,...
        'collectoutput',1);
    cov.(sprintf('%s',types{i1})) = temp{:};
end

fclose(fid);

means_column_info = {'Dust/Total';...
    'SS/Total';...
    'TC/(TC+SO4)';...
    'OC/BC'};

if nvars>4
 means_column_info = {'Dust/Total';...
    'SS/Total';...
    'TC/(TC+SO4)';...
    'OC/BC';'log(extR)'};
end

