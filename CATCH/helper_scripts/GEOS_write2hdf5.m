% This file just cleans up the datafiles to combine them into one .h5 file.
% .h5 is newer and better and will be using that instead. Kind of annoying
% because I have to change everything :-/.
% This sometimes has to be modified accordingly for what you want to do, just FYI.
% One directory holds the SOA tracers, and one holds the other OPs.
clear; clc;

% fp = '/Volumes/cloud/GEOSChem/AEROSOL_Data/Optical_HOURLY/SOA_tracers/';
% fp = '/Volumes/cloud/GEOSChem/AEROSOL_Data/Optical_HOURLY/OP/';
%fp = '/Volumes/cloud/GEOSChem/AEROSOL_Data/Optical_HOURLY/AOD/';
fp = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/geoschem_variables/newemissions/';

Title = 'WARNING!';
STRNG = 'THIS SCRIPT MODIFIES FILE CONTENTS...CONTINUE?';
%choice = questdlg(STRNG,Title);
choice = 'yes';
choice = strcmpi(choice,'yes');

if choice

% I am modifying this to read .mat files now.    
allvars = dir([fp,'*.hdf']);
%allvars = dir([fp,'*.mat']);
allvars = {allvars.name}';
dates = regexp(allvars,'201\d[\d]+','match');
xx = @(x) x{:};
dates = cellfun(xx,dates,'UniformOutput',0);
udates = unique(dates);
nf = numel(udates);

    for j1 = 1:nf
        idx = strcmp(udates(j1),dates);

        h5file = ['Global_DAQ_Flights_HOURLY_0.25x0.3125.bpch.',udates{j1},'.h5'];
        fprintf('Writing to file: %s\n',h5file)
        ish5file = logical(exist(h5file,'file'));
        vars = allvars(idx);

        for i1 = 1:length(vars)
            fn = vars{i1};
            idx = regexp(fn,'[_]*ts','once');
            finfo = hdfinfo([fp,fn]);
            datasets = {finfo.SDS.Name}';
            varname = datasets{end};

            if i1 == 1 && ~ish5file
                for i2 = 1:length(datasets)
                    data = hdfread([fp,fn],datasets{i2});
                    varname = datasets{i2};
                    datainfo = whos('data');
                    h5create(h5file,...
                        ['/',varname],datainfo.size,'Datatype',datainfo.class);
                    h5write(h5file,...
                        ['/',varname], data);
                end
            else
                data = hdfread([fp,fn],datasets{end});
                varname = datasets{end};
                datainfo = whos('data');
                h5create(h5file,...
                    ['/',varname],datainfo.size,'Datatype',datainfo.class);
                h5write(h5file,...
                    ['/',varname], data);

            end
        end
    end

else
    fprintf('File: %s ABORTED.\n',mfilename)
end


% EOF
