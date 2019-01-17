% Batch SPM12 normalization estimation
% 17-02-2017 Jonathan Wirsich / Connectlab
function c2_normalizeest(sess_dir, spm12_dir)

    files = spm_select('FPlist', sess_dir, '^ra.*\.nii');
    filecell = cellstr(files);

    clear matlabbatch
    matlabbatch{1}.spm.spatial.normalise.est.subj.vol = {filecell{1}};
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.tpm = {[spm12_dir '/tpm/TPM.nii']};
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.est.eoptions.samp = 3;
    spm_jobman('run',matlabbatch);

end