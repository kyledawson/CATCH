% Getting the phase function, the lidar ratios, the single scatter alb.
% and the hygroscopic scale factors from the input jv_spec_aod.dat file.
% Author:   Kyle Dawson
clear;close all;clc

cmap = colormap(jet(6));

idx = ones(1,42);
idx(8:8+6)=2;
idx(15:15+6)=3;
idx(22:22+6)=4;
idx(29:29+6)=5;
idx(36:36+6)=6;

%% Set the radius vector:
x = logspace(-3,log10(50),1000); dx = diff(log10(x)); dlogx = dx(1); 

%% Size parameter:
% 
% $$\chi=\frac{2\pi r}{\lambda_{550}}$$
spar = 2*pi*x./.550; % Size parameter at lambda = 550 nm

%% Function to output data from GEOS-Chem jv_spec_aod.dat optics file:
% [headers,data,info] = read_OpticalProps_GeosChem('jv_spec_aod.dat',true);
[headers,data,info] = read_OpticalProps_GeosChem('D:\Armazi\GEOSChem\jv_spec_aod.txt',true); % bs edited

% geometric standard deviation for size distribution of components
sd = info.Size_Distribution;

% NOTE: THIS USES THE EFFECTIVE OPTICAL RADIUS! 
%   This collapses the size distribution to a monomodal one represented
%   by the effective optical radius.
sd(:,1)=data(:,3);
sd(isnan(sd(:,2)),2) = 1.6;

% Index of refraction (IOR):
m = info.N_real-info.N_imaginary*1i; 

% Set the dust IOR:
m(1:7) = 1.58+0.014i; % Martin et al. 2002

% Getting the real component
mr = real(m);

% The cos(theta) integral range i.e. 0:pi
mu=linspace(-1,1,10000);
c170 = cosd(170);

% Interested in bacskcatter so going to integrate the phase function from
% 170:180
intrange = mu<c170;

% Getting chemical components: SO4, OC, BC, Dust, SSaccum, SScoarse
types = info.Type;

% Effective optical radius (um):
reff = data(:,3);

% Extinction efficiency output from GChem Mie code and saved in
% jv-spec.dat.
qeff = data(:,2);
clc

% Setting the onscreen output headers:
fprintf(['%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, ',...
    '%s at lambda = 550 nm\n'],...
    'Index', 'Type','r-eff (um)','Sp-phase(sr)',...
    'Sp-mie(sr)','Sp-particle','SSA-phase', 'SSA (mie)', 'SSA-eff',...
    'N-real', 'Color Ratio','f(RH)','Scale Factor');

% Set a figure for phase function:
set(gca,'fontsize',14,'tickdir','out','ticklength',...
    get(gca,'ticklength')*2.5)
box off
ylabel('p(\tau, cos\theta)');xlabel('Scatter angle, \theta')

% These RH are dry at RH=0;
rhdryidx = 8:7:size(data,1);

%% Reconstructing the phase function
% $$p(r,\cos\theta) = f(\Pi_{1:8})$ where $\Pi_{1:8}$ are the GEOS-Chem
% Legendre coefficients:
for i1 = 1:size(data,1)   
% Grabbing the legendre coefficients:
B = data(i1,5:end)';

% Expanding phase function by calculating legendre polynomial from
% legndre coefficients for each component:
for i2 = 1:numel(B)
    Pl(i2,:) = legendreP(i2-1,mu).*B(i2);
end

% add to get the full phase function:
pcos(i1,:)=sum(Pl,1);

% This is a check:
% Must satisfy the nomrmalization condition so that integral p(cos0) = 2;
% and it does.

%% Lidar ratio method 1: Calculate lidar ratio as 
% $$S_p = \frac{4pi}{SSA \times p(r,\mu=-1)}$$
% 
sp(i1,1) = (4*pi)./(data(i1,4)*mean(pcos(i1,intrange))); % units, (sr)

%% The size distributions
% Create a size distribution. Number is not important as each distribution
% gets normalized to a peak number = 1:
coef = [300,sd(i1,1),sd(i1,2)];

% The only exception for size distributions is dust which has a gamma pdf.
if i1 <= 7
    % 2 is alpha, based on eff variance = 0.2 (martin et al. 2002)
    sizedist(i1,:) = pdf('gam',x,2,sd(i1,1)); 
end

if i1 < 7
    continue
elseif i1 == 7
    % Weighting the size distribution for dust by the mass weights of
    % each size bin.
    sizedist = sizedist./repmat(max(sizedist,[],2),1,size(sizedist,2));
    sizedist = .06*(sum(.25.*sizedist(1:4,:),1))+.12*sizedist(5,:)+...
        .24*sizedist(6,:)+.58*sizedist(7,:);
    
    % This is an alternative and is currently my method of choice.
    % Overwrite the dust size distribution with fit data from Curci:
    fun = @(coef,x) ...
        coef(3)*x.^((1-3*coef(1))/coef(1)).*exp(-x./(coef(1)*coef(2)));
    coef = [0.2477    0.7145   16.3559];
    sizedist = fun(coef,x);
    sizedist = sizedist./max(sizedist);
else
    sizedist = dNdlogDp(coef,x);
end

sizedist = sizedist./max(sizedist); % normalized size dist.

%% Lidar ratio method 2: Calculate as extinction/backscatter coefficient
% $$S_p = \frac{\sigma}{\beta}$$
%
for i2 = 1:length(spar)
    [~,~,qext(i2),qsca(i2),qback(i2),~]=bhmie(spar(i2),m(i1),10);
    
    % This is at 1064 nm, WARNING! IOR is not changed.
    % bhmie is the Bohren and Huffman Mie Code available at link below:
    % http://www.ece.neu.edu/courses/eece5698/2011sp/matlab/bhmie.m
    [~,~,qext2(i2),qsca2(i2),qback2(i2),~]=...
        bhmie(spar(i2).*(550/1064),m(i1),10);
end

% Convert dNdlogDp to dNdDp
sizedistN = sizedist./(2.303.*x);

% Integrating with the size distribution to get extinction and backscatter
% coefficients:
ext = trapz(x,qext.*pi.*x.^2.*sizedistN); % This is per Mm
sca = trapz(x,qsca.*pi.*x.^2.*sizedistN);
bks = trapz(x,qback.*pi.*x.^2.*sizedistN);
bks2 = trapz(x,qback2.*pi.*x.^2.*sizedistN);

if any(i1 == [1:7,rhdryidx])
    extd = ext;
    reff_dry = reff(i1);
    qeff_dry = qeff(i1);
end

% An attempt to get the color ratio. This is not really valid as my index
% of refraction data is only defined at 550 nm. Ignore this output for now.
cr(i1,1) = bks./bks2;

% Method 2: 
sp2(i1,1) = ext./bks;
ssa(i1,1) = sca./ext;

%% Lidar ratio method 3: per particle; extinction to backscatter xsection
cross_idx = closest(spar,sd(i1,1)*2*pi/.550);
sp_cross(i1,1) = qext(cross_idx)./qback(cross_idx);
ssa_cross(i1,1) = qsca(cross_idx)./qext(cross_idx);

% Hygroscopic scale factor:
if i1 <= 7
    fRH(i1,1) = 1;
    scaleF(i1,1) = 1;
else
    fRH(i1,1) = ext./extd;
    scaleF(i1,1) = qeff(i1)./qeff_dry.*(reff(i1)./reff_dry).^2;
end

clear qext qsca qback

% Send output to screen:
fprintf('%d, ',i1);
fprintf('%s, ',types{i1});
fprintf('%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',[data(i1,3),...
    sp(i1),sp2(i1),sp_cross(i1),data(i1,4),ssa(i1),ssa_cross(i1),...
    mr(i1),cr(i1),fRH(i1),scaleF(i1)]);

% Plotting:
% Don't bother plotting if the phase function has unphysical results:
if ~any(pcos(i1,:)<0)
    semilogy(acos(mu)*180/pi,pcos(i1,:),'color',cmap(idx(i1),:))
    hold on
    drawnow
end

clear Pcos

% Format the plot
if i1 == 1
    set(gca,'fontsize',14,'tickdir','out','ticklength',...
        get(gca,'ticklength')*2.5)
    box off
    ylabel('p(\tau, cos\theta)');xlabel('Scattering angle, \theta')
    set(gcf,'color','w')
end

end


