% This file shows how to perform a batch of +/- comparisons

% load the color model
load('../template_colormodel/CM120312.mat');

% set up metadata
experimentName = 'LacI Transfer Curve';
device_name = 'LacI-CAGop';
inducer_name = '100xDox';

% Configure the analysis
% Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
bins = BinSequence(4,0.1,10,'log_bins');

% Designate which channels have which roles
input = channel_named(CM, 'EBFP2');
output = channel_named(CM, 'EYFP');
constitutive = channel_named(CM, 'mKate');
AP = AnalysisParameters(bins,{'input',input; 'output',output; 'constitutive' constitutive});
% Ignore any bins with less than valid count as noise
AP=setMinValidCount(AP,100');
% Ignore any raw fluorescence values less than this threshold as too contaminated by instrument noise
AP=setPemDropThreshold(AP,5');
% Add autofluorescence back in after removing for compensation?
AP=setUseAutoFluorescence(AP,false');
% By default, analysis tries to fit constitutive to transformed and non-transformed components
% If your distribution is more complex or less complex, you can change the number of components
% AP=setNumGaussianComponents(AP,3);

% Make a map of the batches of plus/minus comparisons to test
% This analysis supports two variables: a +/- variable and a "tuning" variable
stem1011 = '../example_assay/LacI-CAGop_';
batch_description = {...
 {'Lows';'BaseDox';
  % First set is the matching "plus" conditions
  {0.1,  {[stem1011 'B9_B09_P3.fcs']}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {[stem1011 'B10_B10_P3.fcs']}};
  % Second set is the matching "minus" conditions 
  {0.1,  {[stem1011 'B3_B03_P3.fcs']};
   0.2,  {[stem1011 'B4_B04_P3.fcs']}}};
 {'Highs';'BaseDox';
  {10,   {[stem1011 'C3_C03_P3.fcs']};
   20,   {[stem1011 'C4_C04_P3.fcs']}};
  {10,   {[stem1011 'B9_B09_P3.fcs']};
   20,   {[stem1011 'B10_B10_P3.fcs']}}};
 };

% Execute the actual analysis
OSbin = OutputSettings('',device_name,'','plots/');
results = process_plusminus_batch( CM, batch_description, AP, OSbin);

% Make additional output plots
for i=1:numel(results)
    OS = OutputSettings(batch_description{i}{1},device_name,'','plots/');
    OS.PlotTickMarks = 1;
    plot_plusminus_comparison(results{i},OS)
end

save('-V7','LacI-CAGop-plus-minus.mat','batch_description','AP','results');
