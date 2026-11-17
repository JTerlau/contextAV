

function plot_networks(env)


[surf_lh, surf_rh] = load_conte69();

labeling = load_parcellation('schaefer',400);

%% MD parcels
% load schaefer atlas
atlas = ft_read_atlas('/Volumes/helfrich_data/ieeg_quenching/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii');


%% just for visualization! this is an approximation with imperfect dice overlap, taken from
%% https://mvpa.blogspot.com/2021/09/approximately-matching-different.html
md_parcels = [99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, ...
                   163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, ...
                   340, 345, 349, 350, 351, 352, 354, 361, 365, 387];

%% FPCN and CON parcels
% Read the label file (tab-separated: index, name, R, G, B, A)
schaefer_path = '/Volumes/jterlau/iEEG_contextA_V/Schaefer2018_400Parcels_17Networks_order.txt';
fid    = fopen(schaefer_path);
labels = textscan(fid, '%d %s %d %d %d %d', 'Delimiter', '\t');
fclose(fid);

parcel_idx   = labels{1};
parcel_names = labels{2};

% Extract parcel numbers by network name string
con_parcels  = parcel_idx(contains(parcel_names, 'SalVentAttn'));
fpcn_parcels = parcel_idx(contains(parcel_names, 'Cont'));


% create masks from Schaefer atlas
mask.md   = double(ismember((1:400)', md_parcels));
mask.fpcn = double(ismember((1:400)', fpcn_parcels));
mask.con  = double(ismember((1:400)', con_parcels));


% replace zeros with NaN
labeling.schaefer_400(labeling.schaefer_400 == 0) = NaN;

networkLabels = {'md','fpcn','con'};

for iNet = 1:3

    thisNet = networkLabels{iNet};

    h = plot_hemispheres(mask.(thisNet), {surf_lh,surf_rh}, ...
        'parcellation',labeling.schaefer_400);


    cm = [0.95 0.95 0.95;[9 74 180]./256];
    colormap(cm)

    % set figure size
    set(gcf,'Position',[631 21 765 775])


    %% save
    cd(fullfile(env.path.main))

    savename = ['brainplot_' thisNet '_network'];
    print(gcf,savename,'-djpeg');
end
