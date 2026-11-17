
function result = electrodes_in_networks(elec_mni,nii_sch,nii_md)

% elec_mni : Nx3 matrix of MNI coordinates [x y z]
% nii_md   : loaded MD atlas (ft_read_mri output)
% nii_sch  : loaded Schaefer atlas (ft_read_mri output)

schaefer_txt = % [path to Schaefer2018_400Parcels_17Networks_order.txt]
dist_thresh  = 5;  % mm
thresh_vox = dist_thresh / 2;  % 5mm / 2mm per voxel = 2.5 voxels

%% --- Build binary masks ---
vol_md  = double(nii_md.anatomy) > 0;
vol_sch = double(nii_sch.anatomy);

% verify both atlases are same resolution and dimensions
assert(isequal(size(vol_md), size(vol_sch)), 'Atlas dimensions do not match')
vox_size = sqrt(sum(nii_md.transform(1:3,1:3).^2));
% fprintf('Voxel size: %.1f x %.1f x %.1f mm\n', vox_size)

% load Schaefer labels
fid          = fopen(schaefer_txt);
labels       = textscan(fid, '%d %s %d %d %d %d', 'Delimiter', '\t');
fclose(fid);
parcel_idx   = double(labels{1});
parcel_names = labels{2};

vol_con  = ismember(vol_sch, parcel_idx(contains(parcel_names, 'SalVentAttn')));
vol_fpcn = ismember(vol_sch, parcel_idx(contains(parcel_names, 'Cont')));

%% --- Precompute non-zero voxel positions ---
[vx,vy,vz] = ind2sub(size(vol_md),   find(vol_md));   md_vox   = [vx,vy,vz];
[vx,vy,vz] = ind2sub(size(vol_con),  find(vol_con));  con_vox  = [vx,vy,vz];
[vx,vy,vz] = ind2sub(size(vol_fpcn), find(vol_fpcn)); fpcn_vox = [vx,vy,vz];

%% --- Transform MNI -> voxel space (single transform) ---
nElec  = size(elec_mni, 1);
elec_h = [elec_mni, ones(nElec, 1)]';
vox    = (inv(nii_md.transform) * elec_h)';  % N x 4
vox    = vox(:, 1:3);                         % N x 3

%% --- Assign electrodes ---
md_dist   = nan(nElec, 1);
con_dist  = nan(nElec, 1);
fpcn_dist = nan(nElec, 1);

for i = 1:nElec
    md_dist(i)   = sqrt(min(sum((md_vox   - vox(i,:)).^2, 2)));
    con_dist(i)  = sqrt(min(sum((con_vox  - vox(i,:)).^2, 2)));
    fpcn_dist(i) = sqrt(min(sum((fpcn_vox - vox(i,:)).^2, 2)));
end

%% --- Output ---
result.md.inside   = md_dist   <= thresh_vox;
result.md.dist     = md_dist;
result.con.inside  = con_dist  <= thresh_vox;
result.con.dist    = con_dist;
result.fpcn.inside = fpcn_dist <= thresh_vox;
result.fpcn.dist   = fpcn_dist;
end





