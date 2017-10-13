function [aid,probability,savepath] = output_aerosol_id(matfilename,training_set_version)
% This script will output aerosol types using the GEOS-Chem means and covariance matrices.
% There are a few notes here:
% 1) There are two outputs. First is aerosol type, second is probability of each class.
% 	The probability of each class is stored in a structure as "probability.<class>"
% 2) The dust probability is calculated as P(total) = P(mahal) * P(dust_model);
% 	This is simply the mahalanobis distance converted to probability and multiplied by
%	dust fraction AOD. This further constrains the dust cluster.
% 3) Mahalanobis distance is converted to chi square cumulative probability. Strictly
% 	speaking this is the probability of the datapoint being random. Therefore,
%	we convert this to probability that the datapoint belongs to the cluster
% 	as follows: P(mahal) = 1-P(random_mahal).
% 4) The output filename is the input filename with an appended "AID" 
%	(i.e. <input_filename>_AID.mat
% 5) All outputs and inputs are in matlab file format. A script to write these to hdf5 
%	format may be used and is inside this package. 
% 6) Training set version is overwritten to use the latest. You can change this.
% Kyle Dawson
%
% UPDATES:
%
% Note that dust probability no longer considers the dust fraction constraint.
% The cluster was built with dust fraction > 0.85 and no outlier trimming. - KD 02-Jan-2017 
% See below for information on the factor flag.
%
% The factor flag allows you to specify an increase in the covariance matrix. If you supply factor = 4,
% then the covariance matrix becomes a 2sigma*2sigma instead of a sigma*sigma. Currently, no 
% factor <fact> flag is used and the covariance matrix is just 1-sigma.

fact = 1;

[mpath,mfile] = fileparts(matfilename);

if strcmp(mpath,'')
	mpath = ['/gpfs_backup/meskhidze_data/DISCOVER_AQ/cluster_variables/newemissions/'];
end

matfilename = [mpath,filesep,mfile,'.mat'];

% Loading the data:
load(matfilename)

% Get training set:
training_set_version = 1;
[means,cov,means_column_info] = read_prespecified_clusters_GC(training_set_version);

fields = fieldnames(cov);
nf     = numel(fields);

% Construct the data matrix for getting MD.

n  = numel(optics.oa2bc);
v1 = reshape(optics.c2cso4,n,1);
v2 = reshape(optics.oa2bc,n,1);
v3 = reshape(optics.sp_eff,n,1);
v4 = reshape(optics.ssa_eff,n,1);
v5 = reshape(optics.nr_eff,n,1);
s = size(optics.oa2bc);
data = [v1,v2,v3,v4,v5];

priordust = optics.dustfraction;
priordustv = reshape(priordust,n,1);
prob 	  = zeros(n,5);

% Order of types is 1, Dust
% 		    2, Marine
%		    3, Urban
%		    4, Smoke
%		    5, Fresh Smoke
for i1 = 1:nf
	 %dbstop in output_aerosol_id at 71
	[MD,P] = mahalanobis(data,fact.*cov.(fields{i1}),means(i1,:));
	
	prob_mahal(:,i1) = 1-P;

	prob_belongs.(fields{i1}) = reshape(100.*(1-P),s);
	
	if false %i1 == 1
		prob_belongs.(fields{i1}) = 100.*(prob_belongs.(fields{i1})./100.*priordust);
	end
end

% find types:
% Add in the dust fraction: This isn't clear to me that it is needed. So I took it out.
prob_mahal(:,1) = prob_mahal(:,1);%.*priordustv;

[maxP,maxPidx] = max(prob_mahal,[],2);

normP = maxP./sum(prob_mahal,2);

rmv = isnan(normP); % This means all probabilities were == 0 and no classification is attempted, outlier or otherwise.

maxP(rmv) = NaN;
maxPidx(rmv) = NaN;

isoutlier = maxP < 0.001;  % i.e. > 99.9 percent chance of being random
isoverlap = normP < 0.6 & ~isoutlier;

aid = maxPidx;
aid(isoutlier) = 6;
aid(isoverlap) = 7;
aid(rmv) = 0;
aid = uint8(aid);
aid = reshape(aid,s);

probability = prob_belongs;

savepath = [mpath,'../../aerosol_types/newemissions/',filesep,mfile,'_aid.mat'];
fprintf('saving to: %s...',savepath)
fprintf('done.\n')






 
