% extract signal from Shirer atlas [Shirer 2012, Cerbral Cortex]
% Import from Sepideh connectivoty scripts
% 17-02-2017 Jonathan Wirsich / Connectlab
function e_extract_ts_shirer(img_dir, atlas_dir)

    img_filter = '^raf.*\.nii$';
    images = cellstr(spm_select('List',img_dir,img_filter));
    %             if isempty(images)
    %                 img_filter = '^w.*\.nii$';
    %                 images = cellstr(spm_select('List',img_dir,img_filter));
    %             end
    %             print(['Found ' num2str(size(images,1)) 'volume files for subject ' subjectList{s}]);
    imgs = char(cellstr(strcat(img_dir,filesep,images)));
    % size(imgs)

    shirer = readShirer([atlas_dir 'shirer' filesep]);
    
    V = spm_vol(imgs);
    timeseries = zeros(2, length(V));
    
   maskfile_gm = cellstr(spm_select('List',img_dir,'^c1.*\.nii$'));
   VGMmask = char(cellstr(strcat(img_dir,filesep,maskfile_gm(1))));
   VGMmask = spm_vol(VGMmask);
    
    count = 0;
    for j = 1:length(shirer)
        display(['Extracting RSN - ' shirer(j).name])
        for it = 1:length(shirer(j).idx)
            count = count+1;
            maskfile = fullfile([img_dir 'shirer' filesep shirer(j).name filesep num2str(it, '%02d') ... 
                filesep 'w' int2str(it) '.nii']);

            maskfile = fullfile(maskfile);
            Vmask = spm_vol(maskfile);

   %extract ts
   %% Function - signal_calc
%----------------------------------------------------
% function sig = signal_calc(V,maskfile,maskname,GMmaskfile,DoGMmask)
% Integrate the values in an segmented image.

    
    fprintf('\n RSN %s - %d: Calculating signal per image...',shirer(j).name, it);
    f = 'img.*(maskI>0);';     
%     f = 'img.*(maskI>0&maskII>0.1);';       
%     if DoGMmask
%         f = 'img.*(maskI>0.5 & maskII>0.5);';       
%     else
%         f = 'img.*(maskI>0.5);';
%     end

    sig = zeros(length(V),1);
    vox = zeros(length(V),1);
    spm_progress_bar('Init',length(sig),'calc signal','images completed');
    fprintf('\n1    ');
    for i = 1:length(sig)      % for each volume (img-file)
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
        timeseries(count, :) = sig;
    end
    spm_progress_bar('Clear');
    fprintf('\n');
    
        end
    end 
    %save timecoursematrix
    save([img_dir filesep 'timeseries_raw.mat'], 'timeseries')
end
