% extract global signal from SPM12 greymatter segmentation on T*2 image
% Import from Sepideh connectivoty scripts
% 17-02-2017 Jonathan Wirsich / Connectlab
function e_extract_ts_globalsig(img_dir)

    img_filter = '^raf.*\.nii$';
    images = cellstr(spm_select('List',img_dir,img_filter));

    imgs = char(cellstr(strcat(img_dir,filesep,images)));
  
    V = spm_vol(imgs);
    timeseries = zeros(1, length(V));
    
   maskfile = cellstr(spm_select('List',img_dir,'^c1.*\.nii$'));
   Vmask = char(cellstr(strcat(img_dir,filesep,maskfile(1))));
   Vmask = spm_vol(Vmask);

   %extract ts
%% Function - signal_calc
%----------------------------------------------------
% function sig = signal_calc(V,maskfile,maskname,GMmaskfile,DoGMmask)
% Integrate the values in an segmented image.
    
    f = 'img.*(maskI>0.3);';       

    sig = zeros(length(V),1);
    intense = cell(length(V),1);
    vox = zeros(length(V),1);
    spm_progress_bar('Init',length(sig),'calc signal','images completed');
    fprintf('\n1    ');
    for i = 1:length(sig)      % for each volume (img-file)
        intense{i} = double.empty;
        for z = 1:V(i).dim(3)  % for each slice
            img   = spm_slice_vol(V(i),...
                spm_matrix([0 0 z]),V(i).dim(1:2),0);
            maskI = spm_slice_vol(Vmask,...
                spm_matrix([0 0 z]),Vmask.dim(1:2),0);
%             maskII = spm_slice_vol(VGMmask,...
%                 spm_matrix([0 0 z]),VGMmask.dim(1:2),0);
            img = eval(f);      % perform masking
            sig(i) = sig(i) + sum(img(:));              % sum of signal values included
            vox(i) = vox(i) + length(find(img(:)>0));   % number of voxels included
            intense{i} = [intense{i}; img(find(img(:)>0))];
        end;
        sig(i) = sig(i)/vox(i);

        spm_progress_bar('Set',i);
        fprintf('.');
        if rem(i,50)==0
            fprintf('\n');
            fprintf('%d ',i);
            if i < 10; fprintf(' '); end;
            if i < 100; fprintf(' '); end;
            if i < 1000; fprintf(' '); end;
        end;
        
        %put into timcourse matric
        timeseries(1, :) = sig;
    end
    spm_progress_bar('Clear');
    fprintf('\n');
    
    %save timecoursematrix
    save([img_dir filesep 'intense'], 'intense')
    save([img_dir filesep 'globalmean'], 'timeseries')
end
