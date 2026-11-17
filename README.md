# fig3 — Network assignment and visualization

Code for assigning iEEG electrodes to canonical cortical networks and testing
whether anticipatory high-frequency activity scales with network membership.

## Scripts

| Script | Purpose |
|---|---|
| `electrodes_in_networks.m` | Assigns each electrode to the MD, CON (SalVentAttn), and FPCN (Cont) networks by nearest-voxel distance in MNI space |
| `test_networks.m` | Loads atlases, calls `electrodes_in_networks`, fits a linear mixed-effects model of anticipatory HFA difference against network membership, and plots results |
| `lookup_reg.m` | Labels each electrode with its Brainnetome and Schaefer parcel for anatomical reference |
| `plot_networks.m` | Renders network masks on the Conte69 surface using BrainSpace |

## Required atlases

### Multiple-Demand (MD) network mask
- **File:** `allParcels-MD-HE197.nii` (whole MD) or `MDnetwork.nii` (core only)
- **Source:** Fedorenko Lab, MIT — [https://www.evlab.mit.edu/resources-all/download-parcels](https://www.evlab.mit.edu/resources-all/download-parcels)
- Derived from fMRI in 197 participants; available in MNI152 2 mm space.

### Schaefer 400-parcel, 17-network atlas
- **Files:**
  - `Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii` (for `lookup_reg.m` and `plot_networks.m`)
  - `Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii` (for `test_networks.m`)
  - `Schaefer2018_400Parcels_17Networks_order.txt` (parcel label table)
- **Source:** Yeo Lab — [https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations)
- The CON mask uses parcels whose names contain `SalVentAttn`; the FPCN mask uses parcels containing `Cont`.

### Brainnetome Atlas (BNA)
- **File:** `BNA_MPM_thr25_1.25mm.nii`
- **Source:** Bundled with FieldTrip (`fieldtrip/template/atlas/brainnetome/`), or directly from [http://atlas.brainnetome.org](http://atlas.brainnetome.org)
- Used in `lookup_reg.m` for anatomical labeling of electrodes only; not used for network assignment.

## Dependencies

- **FieldTrip** (tested with `fieldtrip-20241025`) — `ft_read_mri`, `ft_read_atlas`
  [https://www.fieldtriptoolbox.org/download](https://www.fieldtriptoolbox.org/download)
- **BrainSpace** — `load_conte69`, `load_parcellation`, `plot_hemispheres` (used in `plot_networks.m`)
  [https://brainspace.readthedocs.io](https://brainspace.readthedocs.io)

## Usage

1. Download the MD and Schaefer atlases and note their paths.
2. In `test_networks.m`, set `md_file` and `schaefer_file` to the 2 mm NIfTI files.
3. In `electrodes_in_networks.m`, set `schaefer_txt` to the path of the Schaefer label `.txt` file.
4. Run `test_networks.m`. It calls `electrodes_in_networks` per subject, fits the LME model, and produces the coefficient plot and electrode scatter.

The electrode-to-region lookup (`lookup_reg.m`) runs separately and saves `_reg.mat` files per subject; these are consumed by `test_networks.m`.
