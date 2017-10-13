% Script adds the cf and extinction masks along the curtains via the index variable.
varpath = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/masks/';
d = dir([varpath,'*.mat']);
d2 = dir([varpath,'../curtains/*.mat']);
nf = numel(d);
cfun = @(x) ~isempty(x);

for i1 = 1:nf
fprintf('Working on %d of %d\n',i1,nf)
        load([varpath,d(i1).name])
	sdate = regexp(d(i1).name,'2013\d\d\d\d','match');
        sdate = sdate{:};
        load([varpath,'../clouds/GEOSFP.',sdate,'.A3cld.025x03125.NA.mat']);
        sdatev = regexp({d2(:).name}',sdate,'once');
        sdatev = cellfun(cfun,sdatev);

        idx = find(sdatev);

        %fn = fieldnames(probability);

        if ~isempty(idx)
        for i2 = 1:numel(idx)

                load([varpath,'../curtains/',d2(idx(i2)).name]);

                maskfields = whos('cloud_ext*');

		for xxi = 1:numel(maskfields)
			temp = eval(maskfields(xxi).name);
			temp = temp(index);
			mask.(maskfields(xxi).name) = temp;
		fprintf('still working\n')
		end
	cfrac = clouds(index);

            save([varpath,'../curtains/',d2(idx(i2)).name],'cfrac','mask','-append')

        end
        end

end

