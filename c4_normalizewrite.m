% Batch SPM12 write nromalization
% 17-02-2017 Jonathan Wirsich / Connectlab
function c4_normalizewrite(sess_dir, atlas_dir)

    files = spm_select('FPlist', sess_dir, '^ra.*\.nii');
    filecell = cellstr(files);
    
    %get deformation
    files2 = spm_select('FPlist', sess_dir, '^y_ra.*\.nii');
    filecell2 = cellstr(files2);
    
    %maskfile = fullfile([atlas_dir 'AAL' filesep], 'ROI_MNI_V4.nii');
    
    shirer = readShirer([atlas_dir 'shirer' filesep]);
    mkdir(sess_dir, 'shirer');
    for j = 1:length(shirer)
        display(['Normalizing RSN - ' shirer(j).name])
        mkdir([sess_dir 'shirer'], shirer(j).name)
        for i = 1:length(shirer(j).idx)
            mkdir([sess_dir 'shirer' filesep shirer(j).name], num2str(i, '%02d'))
            
            maskfile = fullfile([atlas_dir 'shirer' filesep shirer(j).name filesep num2str(i, '%02d') ... 
            filesep int2str(i) '.nii']);

            clear matlabbatch
            matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {filecell2{1}};
            matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.vox = [4 4 4];
            matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.bb = [-78 -112 -70
                                                                             78 76 85];
            matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {filecell{1}};
            matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {maskfile};
            matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {[sess_dir 'shirer' filesep ...
                shirer(j).name filesep num2str(i, '%02d')]};
            matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
            matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
            spm_jobman('run',matlabbatch);
        end
    
    end
    
end