% Batch SPM12 timeslice
% 17-02-2017 Jonathan Wirsich / Connectlab
function a_timeslice(TR, nSlices, isInterleaved, sess_dir)

    
    files = spm_select('FPlist', sess_dir, '^f.*\.nii');
    filecell = cellstr(files);

    clear matlabbatch
    matlabbatch{1}.spm.temporal.st.scans = {filecell};
    matlabbatch{1}.spm.temporal.st.nslices = nSlices;
    matlabbatch{1}.spm.temporal.st.tr = TR; 
    %TR = 3.6
    matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nSlices);
    %TA = TR-(TR/nSlices) 3.528

    %interleaved mode
    if isInterleaved == 1
        so = [nSlices:-2:2 (nSlices-1):-2:1];
    else
        so = 1:nSlices;
    end
    matlabbatch{1}.spm.temporal.st.so = so;
    %[50 48 46 44 42 40 38 36 34 32 30 28 26 24 22 20 18 16 14 12 10 8 6 4 2 49 47 45 43 41 39 37 35 33 31 29 27 25 23 21 19 17 15 13 11 9 7 5 3 1];
    matlabbatch{1}.spm.temporal.st.refslice = nSlices/2;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    spm_jobman('run',matlabbatch);
end