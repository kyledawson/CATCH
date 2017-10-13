% Script adds the aid along the curtains via the index variable.
varpath = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/aerosol_types/newemissions/';
d = dir([varpath,'*.mat']);
d2 = dir([varpath,'../../curtains/newemissions/*.mat']);
nf = numel(d);
cfun = @(x) ~isempty(x);

for i1 = 1:nf
fprintf('Working on %d of %d\n',i1,nf)	
	load([varpath,d(i1).name])
	sdate = regexp(d(i1).name,'2013\d\d\d\d','match');
	sdate = sdate{:};
	load([varpath,'../../cluster_variables/newemissions/vars_',sdate,'.mat']);
	df = optics.dustfraction;
	sdatev = regexp({d2(:).name}',sdate,'once');
	sdatev = cellfun(cfun,sdatev);
	
	idx = find(sdatev);
	
	fn = fieldnames(probability);
	
	if ~isempty(idx)
	for i2 = 1:numel(idx)
		
		load([varpath,'../../curtains/newemissions/',d2(idx(i2)).name]); 
		
		aidc = aid(index);
		probc = probability;
		dustfraction = df(index);

		for i3 = 1:numel(fn)
			probc.(fn{i3}) = probc.(fn{i3})(index);
		end

		%gctc = gcAOD(index);
	
		save([varpath,'../../curtains/newemissions/',d2(idx(i2)).name],'dustfraction','aidc','probc','-append')

	end
	end

end

			
			
	
	
