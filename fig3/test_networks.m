



env = run_startup('cluster_local');

md_file       = %[path to allParcels-MD-HE197.nii]
%md_file       = %[path to MDnetwork.nii]
schaefer_file = %[path to Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii]


%% --- Load MD atlas (source: https://www.evlab.mit.edu/resources-all/download-parcels)---
nii_md  = ft_read_mri(md_file);

%% --- Load Schaefer atlas and build binary masks ---
nii_sch  = ft_read_mri(schaefer_file);

%% load elec network masks and minDiff
subvec  = [];
mask    = [];
minDiff = [];
diffA   = [];
diffV   = [];
mni     = [];

for iSub = 1:numel(env.subjects)

    subID = env.subjects{iSub};

    % idx of mod-gen channels
    respIdx = get_chIdx(env,subID);
    chIdx   = respIdx.gen;

    % elecs
    elecs = importdata(fullfile(env.path.data, 'contextV', subID, [subID '_reg.mat']));

    % minDiff
    V = importdata(fullfile(env.path.data,'contextV',subID,[subID '_hfa_diff.mat']));
    A = importdata(fullfile(env.path.data,'contextA',subID,[subID '_hfa_diff.mat']));

    diffA{iSub} = A.diff(chIdx);
    diffV{iSub} = V.diff(chIdx);

    minDiff{iSub} = min([diffA{iSub},diffV{iSub}],[],2);

    % md network binary mask
    result = electrodes_in_networks(elecs.elecpos,nii_sch,nii_md);
    mask.md{iSub}   = result.md.inside(chIdx);
    mask.con{iSub}  = result.con.inside(chIdx);
    mask.fpcn{iSub} = result.fpcn.inside(chIdx);
 
    % subject ID
    subvec = [subvec; repelem(iSub,sum(chIdx),1)];

    % MNI coordinates
    mni{iSub} = elecs.elecpos(chIdx,:);
end

minDiff   = cell2mat(minDiff(:));
diffA     = cell2mat(diffA(:));
diffV     = cell2mat(diffV(:));
mni       = cell2mat(mni(:));
mask.md   = cell2mat(mask.md(:));
mask.con  = cell2mat(mask.con(:));
mask.fpcn = cell2mat(mask.fpcn(:));


%% minDiff vs network masks

tbl = table(double(minDiff),mask.md,mask.fpcn,mask.con,...
    categorical(subvec),'VariableNames',{'minDiff','md','fpcn','con','subject'});

lme = fitlme(tbl,'minDiff ~ md + fpcn + con + (1|subject)');
disp(lme)

%% plot  coefficients
coef = lme.Coefficients;

% extract fixed effects (exclude intercept)
idx    = ~strcmp(coef.Name, '(Intercept)');
names  = coef.Name(idx);
est    = coef.Estimate(idx);
ci_low = coef.Lower(idx);
ci_high= coef.Upper(idx);

% clean labels
names = strrep(names, '_1', '');
names = strrep(names, '_true', '');

% convert CIs to error bar format (distance from estimate)
err_low  = est - ci_low;
err_high = ci_high - est;

figure; hold on
colors = AVcolors;
c = colors.gen.seq;
yline(0, 'k', 'LineWidth', 0.5)
er = errorbar(1:numel(est), est, err_low, err_high,'o');
er.Color = c;
er.MarkerFaceColor = c;
er.MarkerSize = 5;
er.LineWidth = 1;
er.LineStyle = 'none';

xlim([0.5 numel(est)+0.5])
xticks(1:numel(est)); xticklabels(names)
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxisLocation = 'left';
ax.YTickLabelRotation = 0;
set(gcf,'Position',[697 621 326 77])
ylabel('beta')



%% plot scatter

scatter(diffA,diffV,15,minDiff,'filled')
xlabel('Delta HFA aud.')
ylabel('Delta HFA vis.')
axis square
cb = colorbar;
ylabel(cb,'minDiff')
set(gcf,'Position',[792 576 208 122])
c = cbrewer2('YlGnBu');
colormap(c)
lsline

tbl = table(double(diffA),double(diffV),categorical(subvec_adj),'VariableNames',{'diffA','diffV','subject'});

lme = fitlme(tbl,'diffA ~ diffV + (1|subject)');
disp(lme)

%% plot elecs
% make all x coordinates negative to mirror right electrodes to left
% hemisphere

mni_L = mni;
idx_R = mni(:,1) > 0;
mni_L(idx_R,1) = mni(idx_R,1).*-1;

plot_elecs_dfbeta('left','l',mni_L,minDiff)






%% subject averages
uSub = unique(subvec);
nSub = numel(uSub);
subAvg.md   = NaN(nSub,2);
subAvg.con  = NaN(nSub,2);
subAvg.fpcn = NaN(nSub,2);

for iSub = 1:nSub

    subIdx = subvec == uSub(iSub);

    subAvg.md(iSub,1) = mean(minDiff(subIdx & ~mask.md));
    subAvg.md(iSub,2) = mean(minDiff(subIdx & mask.md));

    subAvg.con(iSub,1) = mean(minDiff(subIdx & ~mask.con));
    subAvg.con(iSub,2) = mean(minDiff(subIdx & mask.con));

    subAvg.fpcn(iSub,1) = mean(minDiff(subIdx & ~mask.fpcn));
    subAvg.fpcn(iSub,2) = mean(minDiff(subIdx & mask.fpcn));
end

