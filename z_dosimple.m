%% Do simple Resting state analysis 
% T*2 image only
% 17/02/2016 Jonathan Wirsich / Connectlab
%% configuration
sess_dir = ''; %put folder name containing subject fMRI data; 
atlas_dir = ''; %.../simple-rsfmri/atlas/aal;
spm12_dir = ''; %.../MATLAB/spm12/;

%timeslice
TR = 2.0;
nSlices = 34;
%
a_timeslice(TR, nSlices, 0, sess_dir)
%realign
b_realign(sess_dir)

fprintf('Reorienting ACPC for: %s - %s \n', 'test', sess_dir);
conf.dir = sess_dir;
conf.filepattern = '^af.*\.nii';
reorient_acpc('init', conf);
reply = input('Press Enter to continue', 's');

% marsbar needs to be unloaded because it overwrites spm8 features which
% will be called later
addpath(genpath([spm12_dir 'toolbox' filesep 'marsbar-0.44']));
marsbarTools = SPMmarsbarTools('^af.*\.nii');
%define+extract LCR and white matter
marsbarTools.extractRois(sess_dir);
rmpath(genpath([spm12_dir 'toolbox' filesep 'marsbar-0.44']));

%normalize
c3_segmentnormal(sess_dir, spm12_dir)
c4_normalizewrite(sess_dir, atlas_dir)
% 
% %extract timecourses 
e_extract_ts_shirer(sess_dir, atlas_dir);

e_extract_ts_globalsig(sess_dir);

% %regress
f_regress_nuissance(sess_dir, atlas_dir, 1);

load([sess_dir 'timeseries_regressed.mat']);

filtered = g_filter(regsout, TR);

[reject, accept, param] = scrub('', sess_dir);
scrubbed = filtered(:,accept);

%save connectivity matrix
conn = corr(transpose(scrubbed));
save([sess_dir 'conn_mtrx'], 'conn')

%show scrubbed connectivity matrix
figure; imagesc(conn);
colormap jet

shir = readShirer([sess_dir 'shirer/']);
count = 0.5;
for j = 1:length(shir)
    count = count+length(shir(j).idx);
    shirIdx(j) = count;
    col{j} = 'black';
    label{j} = shir(j).name;
    label2{j} = '';
end
vline(shirIdx, col, label2)
hline(shirIdx, col, label2)
