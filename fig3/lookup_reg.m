

env = run_startup('cluster_local');

%% --- Load atlases ---
bna_path    = %[path to fieldtrip-20241025/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii]
schaefer_path = %[path to Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii];

bna      = ft_read_atlas(bna_path);
schaefer = ft_read_atlas(schaefer_path);

%% --- Precompute non-zero voxel positions for both atlases ---
% Brainnetome
bna_dims  = size(bna.tissue);
[vx, vy, vz] = ind2sub(bna_dims, find(bna.tissue > 0));
bna_vox   = [vx, vy, vz];
bna_vals  = double(bna.tissue(bna.tissue > 0));

% Schaefer
sch_dims  = size(schaefer.parcellation);
[vx, vy, vz] = ind2sub(sch_dims, find(schaefer.parcellation > 0));
sch_vox   = [vx, vy, vz];
sch_vals  = double(schaefer.parcellation(schaefer.parcellation > 0));

%% --- Loop over subjects ---
for iSub = 1:numel(env.subjects)
    subID = env.subjects{iSub}; disp(subID)

    % load electrode coordinates
    elecs = importdata(fullfile(env.path.data,'contextV',subID,[subID '_sa_ranks.mat']));

    % add electrode E2 for TUE13
    if strcmp(subID,'TUE13')
        elecs.eleclabel = [elecs.eleclabel(1:23); 'E2';         elecs.eleclabel(24:end)];
        elecs.elecpos   = [elecs.elecpos(1:23,:); [15.0165	37.2444	32.9853]; elecs.elecpos(24:end,:)];
    end

    nElec    = size(elecs.elecpos, 1);
    elec_mni = elecs.elecpos;

    % --- Transform MNI -> voxel space for both atlases ---
    elec_h = [elec_mni, ones(nElec, 1)]';          

    bna_vox_elec = (bna.transform      \ elec_h)'; % N x 4
    sch_vox_elec = (schaefer.transform \ elec_h)'; % N x 4

    bna_vox_elec = bna_vox_elec(:, 1:3);           % N x 3
    sch_vox_elec = sch_vox_elec(:, 1:3);           % N x 3

    % --- Initialise output ---
    reg.eleclabel       = elecs.eleclabel;
    reg.elecpos         = elec_mni;
    reg.bna_label        = repmat({'unknown'}, nElec, 1);
    reg.bna_parcel       = nan(nElec, 1);
    reg.schaefer_parcel = nan(nElec, 1);
    reg.schaefer_label  = repmat({'unknown'}, nElec, 1);
    reg.bna_all_labels          = bna.tissuelabel';
    reg.schaefer_all_labels     = schaefer.parcellationlabel';

    for iElec = 1:nElec

        % --- Brainnetome: nearest non-zero voxel ---
        dists        = sum((bna_vox - bna_vox_elec(iElec,:)).^2, 2);
        [~, nearest] = min(dists);
        bna_idx      = bna_vals(nearest);
        reg.bna_parcel(iElec) = bna_idx;
        reg.bna_label{iElec}  = bna.tissuelabel{bna_idx};

        % --- Schaefer: nearest non-zero voxel ---
        dists        = sum((sch_vox - sch_vox_elec(iElec,:)).^2, 2);
        [~, nearest] = min(dists);
        sch_idx      = sch_vals(nearest);
        reg.schaefer_parcel(iElec) = sch_idx;
        reg.schaefer_label{iElec}  = schaefer.parcellationlabel{sch_idx};

    end

    % --- Save ---
    for ctx = {'contextA', 'contextV'}
        savepath = fullfile(env.path.data, ctx{1}, subID, [subID '_reg.mat']);
        save(savepath, 'reg')
    end
    clear reg

end



