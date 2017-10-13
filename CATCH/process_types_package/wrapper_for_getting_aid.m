% Wrapper for getting the file ID:
%dbstop if caught error
pth = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/cluster_variables/newemissions/';
d = dir([pth,'*.mat']);
tversion = 1;

for i1 = 1:numel(d)
tstart = tic;	
fprintf('Working on file %d of %d...\n',i1,numel(d))
	[aid,probability,spath] = output_aerosol_id(d(i1).name,tversion);

	save(spath,'aid','probability')
fprintf('Done. ')
telapsed = toc(tstart)./60;
fprintf('(%f minutes elapsed)\n',telapsed);

end
