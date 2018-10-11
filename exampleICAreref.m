
clear;
clc
ft_defaults;

%% loading... 
% load eeg and the samplinginfo field
load('eeg.mat');

cfg = [];
cfg.latency = [-1 1.5];

eeg = ft_selectdata(cfg,eeg);

%% do the ica with a highpass of 1Hz
cfg             = [];
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 1;
tmp             = ft_preprocessing(cfg, eeg);
data_comp       = ft_componentanalysis([], tmp); %using eeglab runica

unmixing        = data_comp.unmixing;
used_labels     = tmp.label;
mixing          = inv(unmixing);

save unmixing unmixing
% mixing * components = data
% unmixing * data = components 

cfg          = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg,data_comp);


% maybe a good point to save the unmixing matrix and the used_labels
% check the distribution of hidden sources;
figure;
imagesc(abs(mixing))
caxis([-max(abs(mixing(:))) max(abs(mixing(:)))]);
colormap(jet(256)); 
colorbar;

set(gca, 'XTick',1:numel(used_labels));
set(gca, 'YTick',1:numel(used_labels));
set(gca, 'YTickLabel', used_labels);
set(gcf, 'Color', 'w');
%% now apply the unmixing matrix to he original eeg (no Highpass filter)
% this is a good point to go on after loading the unmixing matrix
cfg             = [];
cfg.topolabel   = used_labels;
cfg.unmixing    = unmixing;
data_comp       = ft_componentanalysis(cfg, eeg);
%% do a chi-square test on the columns of the unsorted mixing matrix
chi_val = sum(bsxfun(@rdivide, ((bsxfun(@minus,abs(mixing),  mean(abs(mixing),1))).^2), ...
    mean(abs(mixing),1)),1);

%sort the chi-square values and keep track of the sorting-indices
[chisort, s_indices] = sort(chi_val, 'descend');

% now find cutoff values for different critical p-values
crit_01 = numel(find(chisort>chi2inv(0.99, numel(used_labels)-1)));
crit_05 = numel(find(chisort>chi2inv(0.95, numel(used_labels)-1)));
crit_10 = numel(find(chisort>chi2inv(0.9, numel(used_labels)-1)));
crit_20 = numel(find(chisort>chi2inv(0.8, numel(used_labels)-1)));

% define the cutoff for which components are broad
crit_value = crit_20;
%% plot again the mixing matrix(topography) with components sorted by broadness
figure
imagesc((mixing(:, s_indices)))
caxis([-max(abs(mixing(:))) max(abs(mixing(:)))]); %axis xy
colormap(jet(256)); 
colorbar;

set(gca, 'XTick',1:numel(used_labels))
set(gca, 'YTick',1:numel(used_labels))

set(gca, 'YTickLabel', used_labels)
tmp = [1:numel(used_labels)];
set(gca, 'XTickLabel', num2cell(tmp(s_indices)))
colormap(jet(256));

line([crit_value+0.5, crit_value+0.5], [0, numel(used_labels)+0.5] , 'LineStyle', '--', 'LineWidth', 2, 'Color', [1 0 1]);
set(gcf, 'Color', 'w');
%% do bipolar referencing 
 
% create the new labels here
tmplabel        = circshift(eeg.label, -1);
tmp2            = cellfun(@(x, y) strcat(x, ' - ' , y), eeg.label, tmplabel, 'Un', 0);
data_bip.label  = tmp2(1:end-1);

% same time axis etc
data_bip.time       = eeg.time;
data_bip.fsample    = eeg.fsample;
data_bip.trialinfo  = eeg.trialinfo;
data_bip.sampleinfo = eeg.sampleinfo;

% create the matrix that does the bipolar referencing
bip_mat = zeros(numel(data_bip.label), numel(eeg.label));
bip_mat(sub2ind(size(bip_mat), [1:numel(data_bip.label)],...
    [1:numel(data_bip.label)])) = 1;
bip_mat(sub2ind(size(bip_mat), [1:numel(data_bip.label)],...
    [2:numel(data_bip.label)+1])) = -1;
% to check imagesc(bip_mat)

 
% now multply the bip_mat with every trial using cellfun (uniform output
% false)
data_bip.trial = cellfun(@(x) bip_mat*x, eeg.trial, 'Un', 0);
% note: this contains too many channels still. all unplausible channels
% (difference between shafts) need to be rejected manually! see below
% % !!! only for patient1
% sel = true(1,numel(data_bip.label));
% return
% sel([6 12]) = 0;
% cfg         = [];
% cfg.channel = data_bip.label(sel);
% data_bip    = ft_selectdata(cfg, data_bip);

%% do method 4a: reject broad components and backproject
reject_components   = find(~ismember(1:size(data_comp.label), s_indices(1:crit_value)));
cfg                 = []; 
cfg.component       = reject_components;
data_4a             = ft_rejectcomponent(cfg, data_comp);

% cfg = [];
% cfg.viewmode = 'vertical';
% ft_databrowser(cfg,data_4a);

%% do method 4b: keep working with the components
%(i.e. we change the labels to the elec-label where the source is strongest)

% these are the indices of components we want to keep;
sel_comps       = s_indices(1:crit_value);
channel_labels  = repmat({''}, numel(sel_comps), 1);

% loop through components to generate the labels
for sel = 1 : numel(sel_comps)  
    
    % get the channel-idx where the component peaks
    [~, tmp_ind] =  max(abs(mixing(:,sel_comps(sel))));
    % create the temporary label 
    tmp_label = eeg.label{tmp_ind};    
    % if we already have a source on this channel, count up (a,b,c...) 
    tmp2 = tmp_label; counter_ = 0;
    while ismember(tmp_label, channel_labels)
        tmp_label = [tmp2 char('a' + counter_)];
        counter_ = counter_+1;
    end
    % store the label here
    channel_labels{sel} =  tmp_label;
end
datatmp = data_comp;
%overwrite the labels of the selected components
datatmp.label(sel_comps) = channel_labels;

% now select the local components only
cfg         = [];
cfg.channel = channel_labels; % the local components we want to keep
data_4b     = ft_selectdata(cfg, datatmp);

% cfg = [];
% cfg.viewmode = 'vertical';
% ft_databrowser(cfg,data_4b);
%% method 4c: selective remix

% first select only those columns from the mixing matrix that are local
% (keep the ascending order)
mixing_remain = mixing(:, sort(s_indices(1:crit_value), 'ascend'));

% now find wich channel in each column marks the peak (only for the local
% ones) 
[~, peak_indices] = max(abs(mixing_remain),[], 1);

% inds_tmp are the indices in the mixing_remain that mark the values at peak
% e.g. sub2ind(size(mixing_remain), 13,1) is the 13th element in the matrix
inds_tmp = sub2ind(size(mixing_remain), peak_indices, 1:size(mixing_remain,2));

% now write only the peak values in a same-size-matrix of zeros 
new_mix = zeros(size(mixing_remain));
new_mix(inds_tmp) = mixing_remain(inds_tmp);

% ... and reconstruct the final n*n mixing matrix with the correct entries
mix_fin = zeros(size(mixing));
mix_fin(:, sort(s_indices(1:crit_value), 'ascend')) = new_mix;

%% just some plots to view everything went fine: 
figure;
imagesc((mixing(:, s_indices)),[-max(abs(mixing(:))) max(abs(mixing(:)))]); %axis xy
colormap(jet(256)); colorbar;

set(gca, 'XTick',1:numel(used_labels))
set(gca, 'YTick',1:numel(used_labels))

set(gca, 'YTickLabel', used_labels)
tmp = [1:numel(used_labels)];
set(gca, 'XTickLabel', num2cell(tmp(s_indices)))
colormap(jet(256));

line([crit_value+0.5, crit_value+0.5], [0, numel(used_labels)+0.5] , 'LineStyle', '--', 'LineWidth', 2, 'Color', [1 0 1]);
set(gcf, 'Color', 'w');
title('sorted mixing matrix')
figure; 
imagesc((mix_fin(:, s_indices)),[-max(abs(mix_fin(:))) max(abs(mix_fin(:)))]); %axis xy
colormap(jet(256)); colorbar;

set(gca, 'XTick',1:numel(used_labels))
set(gca, 'YTick',1:numel(used_labels))

set(gca, 'YTickLabel', used_labels)
tmp = [1:numel(used_labels)];
set(gca, 'XTickLabel', num2cell(tmp(s_indices)))
colormap(jet(256));

line([crit_value+0.5, crit_value+0.5], [0, numel(used_labels)+0.5] , 'LineStyle', '--', 'LineWidth', 2, 'Color', [1 0 1]);
set(gcf, 'Color', 'w');
title('sorted reduced mixing matrix');

figure; 
imagesc((mix_fin(any(mix_fin~=0,2), :)))
caxis([-max(abs(mix_fin(:))) max(abs(mix_fin(:)))]); %axis xy
colormap(jet(256)); 
colorbar;

set(gca, 'XTick',1:numel(used_labels))
set(gca, 'YTick',1:crit_value)

set(gca, 'YTickLabel', (eeg.label(any(mix_fin~=0,2))))
tmp = [1:numel(used_labels)];
set(gca, 'XTickLabel', num2cell(tmp))
colormap(jet(256));
set(gcf, 'Color', 'w');

title('to be multiplied mixing matrix')
%% now create the eeg 
data_4c = eeg;
% only channels that have non-zero weights somewhere will remain
data_4c.label = eeg.label(any(mix_fin~=0,2));

% now multiply the reduced mixing_mat with every component trial using cellfun
%emty channels are excluded here by not using empty rows of the mix_fin
data_4c.trial = cellfun(@(x) mix_fin(any(mix_fin~=0,2),:)*x, data_comp.trial, 'Un', 0);
 
%% just a quick check in the databrowser
cfg             = []; 
cfg.keeptrials  = 'no';
tl4a = ft_timelockanalysis(cfg, data_4a);
tl4b = ft_timelockanalysis(cfg, data_4b);
tl4c = ft_timelockanalysis(cfg, data_4c);

tlbip = ft_timelockanalysis(cfg, data_bip);

cfg          = [];
cfg.viewmode = 'vertical'; 
ft_databrowser(cfg, tl4c);
%% save different versions

eeg = data_bip;
save eeg_ar_bipolar eeg

eeg = data_4a;
save eeg_ar_ICA_a eeg

eeg = data_4b;
save eeg_ar_ICA_b eeg

eeg = data_4c;
save eeg_ar_ICA_c eeg

%%
