%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries: set up TASBE analytics package:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath('~/cur_proj/SynBioIRAD/TASBE_Analytics/');
addpath('your-path-to-analytics');
% turn off sanitized filename warnings:
warning('off','TASBE:SanitizeName');

colordata = '../../colortest/';
dosedata = '../../plasmidtest/';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples of flow data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pure scatter
fcs_scatter([dosedata 'LacI-CAGop_C4_C04_P3.fcs'],'PE-Tx-Red-YG-A','Pacific Blue-A',0,[0 0; 6 6],1);
fcs_scatter([colordata '07-29-11_EYFP_P3.fcs'],'FITC-A','Pacific Blue-A',0,[0 0; 6 6],1);
% smoothed density plot
data1 = fcs_scatter([dosedata 'LacI-CAGop_C4_C04_P3.fcs'],'PE-Tx-Red-YG-A','Pacific Blue-A',1,[0 0; 6 6],1);
data2 = fcs_scatter([colordata '07-29-11_EYFP_P3.fcs'],'FITC-A','Pacific Blue-A',1,[0 0; 6 6],1);

% Things to notice:
% - look at the size of data1 and data2: there's a *LOT* of points in these samples
% - because there is so much of it, pure scatter graphs are not sufficient for interpreting the data
% - the axes are logarithmic, and variation is evenly distributed on the log scale
% - the data runs up against the axes: there are values less than zero not shown
%   less than zero values come from sensor error
% - low values are quantized, but not round
% - the very highest values saturate, at around 10^5.5 in these files
% - the populations of cells are complex, multimodal, and range widely in observed fluorescence


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What's in an FCS 3.0 file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rawdata hdr data] = fca_readfcs([colordata '07-29-11_EYFP_P3.fcs']);

hdr

% Things to notice:
% - Date and time stamp
% - 50K data points gathered in one minute
% - "par" contains all the information about channels
% - "timestep" says what the units of the time field are
% - Most other information is not very useful

{hdr.par(:).name}

% Things to notice:
% - High-dimensional data!
% - (H)eight, (W)idth, and (A)rea channels
% - FSC, SSC: Forward scatter and side scatter
% - Time

hdr.par(7)

% Things to notice:
% - saturation at 2^18
% - configurability of reported resolution: "decade", "log", "logzero"
%   This can cause real problems!
%   c.f. read_filtered_au early data exclusion

sum(data(:,7)<=0)
% How much FITC is less than zero?  Quite a bit!

data(data(:,7)>260000,7)
% ans = 262143
% How high was this really?  We cannot know...




channels = {};
channels{1} = Channel('', 'Pacific Blue-A', '',234,456);
channels{2} = Channel('', 'PE-Tx-Red-YG-A', '',890,123);
channels{3} = Channel('', 'FITC-A', '',789,345);
CM = ColorModel('','',channels,{{},{},{}},{});

filtered = read_filtered_au(CM,[colordata '07-29-11_EYFP_P3.fcs']);
% Notice that 'filtered' is smaller than 'data': we've dropped the first 25 "units" of time
% It should be in tenths of seconds, but in fact it's the uninterpreted timestep field

CM = set_dequantization(CM,true);
[dequantized hdr] = read_filtered_au(CM,[dosedata 'LacI-CAGop_C4_C04_P3.fcs']);
xc = dequantized(:,10); yc = dequantized(:,11);
pos = xc>0 & yc>0;
figure; smoothhist2D(log10([xc(pos) yc(pos)]),10,[200, 200],[],'image',[0 0; 6 6]);
% quantization gets smoothed out by introduction of small noise
