function optics = get_optics(filename)
% INPUTS: filename - a .h5 file with at minimum the following variables:
%                          dust, oc, ssa, ssc, bc and so4 AOD
%                          rh in %
%                          box height
%                          temperature (optional)
%                          pressure (optional)
% OUTPUTS:  optics - a structure array with the following:
%                    effective lidar ratio
%                    effective single scatter albedo
%                    effective index of refraction (real component)
%                    organic aerosol to black carbon (mass ratio)
%                    total carbon to sum of total carbon and sulfate (mass
%                       ratio)
%                    center altitudes (km)
%                    extinction ratio (total to molecular)
% AUTHOR: Kyle Dawson
%
% TO DO:
% From the GEOS-Chem output, create the following variables:
% organic aerosol to black carbon mass ratio
% total carbon to sum of total carbon and sulfate mass ratio
% effective lidar ratio
% effective single scatter albedo
% effective index of refraction
% center altitudes from 0 meters ASL
% dust fractional optical depth
%
% ASSUMPTIONS
% These variables are calculated under the following assumptions:
% 1) if it can get wet, it will get wet. 
%       Hydrophobic fraction at RH > 0 = 0.
% 2) the dust size distribution is constant shape and does not change.
% 3) I can weight the optical properties by their component optical depths.
%
% NEEDED IMPROVEMENTS:
% I really want to compare this to FlexAOD. 
% FlexAOD will simulate optics for the masses including the hydrophobic
% fractions and will also explicitly consider the dust size bins.

%% Begin File
% Initialize a few checks:
[~,~,chk] = fileparts(filename);
musthave = {'BXHEIGHT','OPBC','OPOC','OPSO4','OPSSa','OPSSc','OPD','RH'};
shouldhave = {'TMPU','PSURF'};
cfun = @(x) ~isempty(x);
approx_molecular = false;
if ~strcmp(chk,'.h5')
    error('Wrong file type.')
else
    info = h5info(filename);
    vars = {info.Datasets.Name}';
    
    x = zeros(size(musthave));
    for i1 = 1:numel(musthave)
        xx = regexp(vars,musthave(i1),'once');
        x(i1) = any(cellfun(cfun,xx));
    end
    
    if ~all(x)
        fprintf('I think you are missing: %s\n',musthave{~x})
        error('Missing a variable, cannot continue')
    end
    
    x = zeros(size(shouldhave));
    for i1 = 1:numel(shouldhave)
        xx = regexp(vars,shouldhave(i1),'once');
        x(i1) = any(cellfun(cfun,xx));
    end
    
    if ~all(x)
        fprintf('You are missing: %s\n',shouldhave{~x})
        fprintf('Warning: an approximation will be used.\n')
        approx_molecular = true;
    end      
end

% Checks are passed.

%% Load the data:
for i1 = 1:numel(musthave)
    
    idx = regexp(vars,musthave(i1),'once');
    idx = cellfun(cfun,idx);
    fprintf('Loading %s...',musthave{i1})
    eval(sprintf('%s=h5read(''%s'',''%s'');',musthave{i1},filename,...
        ['/',vars{idx}]))
    fprintf('Done.\n')
end

if ~approx_molecular
for i1 = 1:numel(shouldhave)
    idx = regexp(vars,shouldhave(i1),'once');
    idx = cellfun(cfun,idx);
    
    fprintf('Loading %s...',shouldhave{i1})
    eval(sprintf('%s=h5read(''%s'',''%s'');',shouldhave{i1},filename,...
        ['/',vars{idx}]))
    fprintf('Done.\n')
end
end

%% Load the helper functions and if needed, approximate molecular data.
% Get the current working directory.
currpath = pwd;
% First get molecular data.
% load([currpath,'GEOSChem/molecular_approximate_profiles.mat']);
load molecular_approximate_profiles.mat
if approx_molecular
    t = molext.T; p = molext.P;
    t = [t,nan]; p = [p,nan];
    TMPU = repmat(t,24,1,202,225);
    PSURF = repmat(p,24,1,202,225);
end

Cs = 3.267e-6; % Hostetler et al., CALIOP ATBD molecular constant [K/hPa/m]
% The total to molecular extinction ratio:
B = (OPBC+OPOC+OPSSa+OPSSc+OPD+OPSO4).*TMPU./...
    (BXHEIGHT.*Cs.*PSURF);

%% GET CENTER ALITUDES:
% The center altitudes:
surf = zeros(24,1,202,225);
z_edges = cat(2,surf,cumsum(BXHEIGHT,2));
k = 1:size(z_edges,2)-1;
z_ctrs = z_edges(:,k,:,:)+0.5*BXHEIGHT;

%% GET OPTICS:
% Load the optical data:
% load([currpath,'GEOSChem/ssa_optics_jvspec_v2.mat']);
load ssa_optics_jvspec_v2.mat
for opt = 1:6
    [IDX.(sprintf('t%d',opt)),lr.(sprintf('t%d',opt))] = mapRH_4D(RH,rhlut,LUT_LR(opt,:));
    ssatemp = LUT_SSA(opt,:); ssa.(sprintf('t%d',opt)) = ssatemp(IDX.(sprintf('t%d',opt)));
    nrtemp = LUT_RRI(opt,:); nr.(sprintf('t%d',opt)) = nrtemp(IDX.(sprintf('t%d',opt)));
    alphatemp = LUT_MEXT(opt,:); alpha.(sprintf('t%d',opt)) = alphatemp(IDX.(sprintf('t%d',opt)));
end

% Grab the fractional Optical Depths
tot = OPBC+OPOC+OPSSa+OPSSc+OPD+OPSO4;
x_dst = OPD./tot;
x_ssa = OPSSa./tot;
x_ssc = OPSSc./tot;
x_so4 = OPSO4./tot;
x_oc = OPOC./tot;
x_bc = OPBC./tot;

% Calculate the optical variables:
ssa_eff = x_dst.*ssa.t1 + x_so4.*ssa.t2 + x_bc.*ssa.t3 + x_oc.*ssa.t4 ...
    + x_ssa.*ssa.t5 + x_ssc.*ssa.t6;

% Index of refraction is tricky:
% It must be weighted by total volume in the bin, not the total AOD, but
% they are related by the mass extinction coefficient.
% The assumption for dust is still used. It is weighted by its size
% distribution extinction to calculate average alpha and ssa.
dens.t1 = 2.6410;
dens.t2 = 1.7;
dens.t3 = 1; 
dens.t4 = 1.8;
dens.t5 = 2.2;
dens.t6 = 2.2;

vd = OPD./(alpha.t1.*dens.t1);
vso4 = OPSO4./(alpha.t2.*dens.t2);
vbc = OPBC./(alpha.t3.*dens.t3);
voc = OPOC./(alpha.t4.*dens.t4);
vssa = OPSSa./(alpha.t5.*dens.t5);
vssc = OPSSc./(alpha.t6.*dens.t6);
vtot = vd+vso4+vbc+voc+vssa+vssc;

nr_eff = (vd.*nr.t1 + vso4.*nr.t2 + vbc.*nr.t3 + voc.*nr.t4 ...
    + vssa.*nr.t5 + vssc.*nr.t6) ./ vtot;

sp_eff = tot ./ (OPD./lr.t1 + OPSO4./lr.t2 + OPBC./lr.t3 + OPOC./lr.t4 ...
    + OPSSa./lr.t5 + OPSSc./lr.t6);

% The following is an addition and is from Sugimoto and Lee 2006
% depolarization ratio is parameterized using dust backscatter mixing ratio (X)
X = (OPD./lr.t1)./(OPD./lr.t1 + OPSO4./lr.t2 + OPBC./lr.t3 + OPOC./lr.t4 ...
    + OPSSa./lr.t5 + OPSSc./lr.t6);
dep = @(X) ((1+0.35)./(X.*0.35)-1).^-1;
dep_eff = dep(X);

%% GET MASS RATIOS:
% The masses are calculated as g/m2
oamass = OPOC./alpha.t4;   % Note I say oa because it is OCp + OCsoa
bcmass = OPBC./alpha.t3;
so4mass = OPSO4./alpha.t2;
% Note that various papers indicate OA:OC ratios
% Some common ones are 2.1 by the Kim et al. paper
% Another is 2.24 by the GEOSChem wiki somewhere
% That is that if you want primary OC, then mOC = tauOC/(alphaOC*2.1);
oa2bc = oamass./bcmass;

c2csulfate = (oamass+bcmass)./(oamass+bcmass+so4mass);

%% Assign the structure:

optics.oa2bc = oa2bc;
optics.c2cso4 = c2csulfate;
optics.sp_eff = sp_eff;
optics.nr_eff = nr_eff;
optics.ssa_eff = ssa_eff;
optics.z_centers = z_ctrs;
optics.z_edges = z_edges;
optics.extratio = B;
optics.approx_molec_flag = approx_molecular;
optics.dustfraction = x_dst;
optics.dep_eff = dep_eff;
% EOF

%% CALLBACK FUNCTIONS
% Map the variables to their RH
function [IDX,NEWVAR] = mapRH_4D(inputRH,rhvec,var)
% INPUTS:       inputRH - MxN array
%               rhvec - 1xP vector
%               var - 1xP vector
% OUT:          idx, newvar, MxN mapped indices and variables.
% enter a matrix of inputRH to get out the mapped index for the LUT
% you can optionally enter a variable same size as the rhvector to be
% mapped to the size of the inputRH.

s1 = size(inputRH,1);
IDX = zeros(size(inputRH));
NEWVAR = IDX;

A = rhvec;
inputRH(inputRH>=100)=99.9999;

for i1 = 1:s1

B = squeeze(inputRH(i1,:,:,:));
%dbstop in mapRH.m at 23 
[A,idx1] = sort(A);
[~,idx2] = histc(B,A);
idx = idx1(idx2);

if nargin>2
newvar = var(idx2);
else newvar = [];
end

IDX(i1,:,:,:)=idx;
NEWVAR(i1,:,:,:)=newvar;

end



    




    





