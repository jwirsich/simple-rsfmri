% Batch SPM12 nromalization
% 17-02-2017 Jonathan Wirsich / Connectlab
function c_normalize(sess_dir)

    files = spm_select('FPlist', sess_dir, '^ra.*\.nii');
    filecell = cellstr(files);

    clear matlabbatch
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {filecell{1}};
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = filecell;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/home/jwirsich/Documents/MATLAB/spm12/tpm/TPM.nii'};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    spm_jobman('run',matlabbatch);

end