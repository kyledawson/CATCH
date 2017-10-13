% This script reads in the clouds and extinction ratio
% to generate two aerosol masks.
% First mask is cloud_ext_mask_th_fraction
%	This is a cloud mask with > threshold cloud fraction
%	and an extinction ratio < threshold;
% Second mask is same as first but a more relaxed
%	> 0.1 cloud fraction. It is called
%	cloud_ext_mask_th0p1_fraction
% Saved as variables in matfile: vars_<filedate>_masks.mat
% I recommend applying this mask to the outlier aerosol id only.
% but that must be done in post-processing.
% Kyle Dawson
mainpath = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/';
d = dir([mainpath,'cluster_variables/*.mat']);

nf = numel(d);

th_e = 1.1;   % Extinction ratio threshold value
th_c = 0.001; % Cloud fraction threshold value

for i1 = 1:nf
	fprintf('Working on %d of %d...',i1,nf)	
	[~,shrt,~] = fileparts(d(i1).name);
	load([mainpath,'cluster_variables/',d(i1).name],'optics');
	load([mainpath,'clouds/GEOSFP.',shrt(6:end),'.A3cld.025x03125.NA.mat']);

	extr = optics.extratio;
	cf = clouds;

	cloud_ext_mask_th_fraction = extr<th_e|cf>th_c;
	cloud_ext_mask_th0p1_fraction = extr<th_e|cf>th_c;

	save([mainpath,'masks/masks_',shrt(6:end),'.mat'],'cloud_ext_*')
	fprintf('Done.\n')
end
clear;clc;exit
