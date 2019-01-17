classdef SPMfslTools
    %SPMFSLTOOLBOX Class to call fsl out of Matlab using system()
    %
    %Version 2-12-2013 CEMEREM Jonathan Wirsich
    %tested for Ubuntu 12, 13 + FSL 5.0
    %forked from multimodal private 18-10-2017 Jonathan Wirsich 
    %ConnectLab
    
properties
    filepattern
    c
    type
end    

methods
        
function obj = SPMfslTools(conf, fp, outtype)
    obj.c = conf;
    
    obj.filepattern = fp;
    obj.type = outtype;
    % set FSL environment
    setenv('FSLDIR','/usr/share/fsl/5.0');  % this to tell where FSL folder is
    setenv('FSLOUTPUTTYPE', 'NIFTI'); % this to tell what the output type would be
end

%% Call FSL/Flirt to register T1 to AAL 
% @param structural subjectpath
% writes the tranfomed T1 in AAL (T12AAL.nii)
% and the applied transformation matrix (T12AAL.mat)
function linearRegristration(obj, subjectpath)
    obj.linearRegristrationGeneric(subjectpath, 'T12AAL', [obj.c.aal_dir 'flip_single_subj_T1.nii'], 0)
end  

%% Call FSL/Flirt to register two comparable images
% @param structural subjectpath
% @param tranformation name
% @param reference image
% @param flag if to switch reference with structural image
% writes the tranfomed T1 in AAL ($transfo .nii)
% and the applied transformation matrix ($transfo .mat)
function linearRegristrationGeneric(obj, subjectpath, transfo, ref, invertReference, outpath)
    
   if nargin < 6
       outpath = subjectpath;
   end
    
   files = obj.getFilesIgnoreCase(subjectpath, obj.filepattern);
   
   %switch reference with target
   if(invertReference==1)
       temp = ref;
       ref = files;
       files = temp;
   end
    
    %build up fsl params
    cmd = 'flirt';

    %ref = [obj.c.aal_dir 'flip_single_subj_T1.nii'];
    cmd = [cmd ' -ref ' ref];
    cmd = [cmd ' -in ' files];
    cmd = [cmd ' -out ' outpath transfo obj.type];
    cmd = [cmd ' -omat ' outpath transfo '.mat'];

    cmd = [cmd ' -dof 12']; %affine 12 params
    %incorrecttly oriented
    cmd = [cmd ' -searchrx -180 180'];
    cmd = [cmd ' -searchry -180 180'];
    cmd = [cmd ' -searchrz -180 180'];

    obj.call_fsl(cmd);
end 

%% inverts transformation of coregistered T1 and AAL
% @param structural subjectpath
% needs T12AAL.mat transformation matrix in folder
% writes AAL2T1.mat
function invertTransform(obj, subjectpath)
    obj.invertTransformGeneric(subjectpath, 'T12AAL', 'AAL2T1');
end

%% inverts transformation of corgesitered T1 and AAL
% @param structural subjectpath
% needs $transfo .mat transformation matrix in folder
% writes $transfoinv .mat
function invertTransformGeneric(obj, subjectpath, transfo, transfoinv)
    cmd = 'convert_xfm';

    cmd = [cmd ' -inverse ' subjectpath filesep transfo '.mat'];
    cmd = [cmd ' -omat ' subjectpath filesep transfoinv '.mat'];
    
    obj.call_fsl(cmd);
end

%% applies transformation from AAL2T1 on AAL atlas
% @param structural subjectpath
% needs AAL2T1.mat,  ROI_MNI_V4.nii AAL atlas file
% writes ROI_AAL2T1.nii with transformed AAL atlas ROIs
function applyTransformation(obj, subjectpath)
    obj.applyTransformationGeneric(subjectpath, [obj.c.aal_dir 'ROI_MNI_V4.nii'], 'AAL2T1');
end

%% applies a calculated transformation on parcelation
% @param structural subjectpath
% needs $transfo .mat,  rois.nii parcelation file
% writes ROI_£transfo .nii with transformed parcelations
function applyTransformationGeneric(obj, subjectpath, pattern, rois, transfo)
    cmd = 'flirt';
    
    files = obj.getFilesIgnoreCase(subjectpath, pattern);

    cmd = [cmd ' -applyxfm'];
    cmd = [cmd ' -init ' subjectpath filesep transfo '.mat'];
    cmd = [cmd ' -ref ' files];
    cmd = [cmd ' -in ' rois];
    cmd = [cmd ' -out ' subjectpath filesep 'ROI_' transfo '.nii'];
    cmd = [cmd ' -interp nearestneighbour'];
    
    obj.call_fsl(cmd);
end

%% applies a calculated transformation on other files
% @param structural subjectpath
% needs $transfo .mat,  $sourcename .nii of foiles to be tranfered
% writes r1$sourcename .nii with transformed parcelations
function applyTransformationToOther(obj, subjectpath, pattern, sourcename, transfo, outpath)
    
    if nargin < 6 
        outpath = subjectpath;
    end
    
    cmd = 'flirt';
    
    files = obj.getFilesIgnoreCase(outpath, pattern);
    sourcefiles = obj.getFilesIgnoreCase(subjectpath, ['^' sourcename '\.nii']);

    cmd = [cmd ' -applyxfm'];
    cmd = [cmd ' -init ' outpath filesep transfo '.mat'];
    cmd = [cmd ' -ref ' files];
    cmd = [cmd ' -in ' sourcefiles];
    cmd = [cmd ' -out ' outpath filesep 'r' sourcename '.nii'];
    cmd = [cmd ' -interp nearestneighbour'];
    
    obj.call_fsl(cmd);
end

%% applies a calculated transformation on other files
% @param structural subjectpath
% needs $transfo .mat,  $sourcename .nii of foiles to be tranfered
% writes r1$sourcename .nii with transformed parcelations
function applyTransformationToOtherGeneric(obj, pathref, ref, inPath, inFile, transfopath, transfofile)
%     outformat = '.nii.gz';
    
    cmd = 'flirt';

    cmd = [cmd ' -applyxfm'];
    cmd = [cmd ' -init ' transfopath transfofile '.mat'];
    cmd = [cmd ' -ref ' pathref ref obj.type];
    cmd = [cmd ' -in ' inPath inFile obj.type];
    cmd = [cmd ' -out ' inPath inFile '_' transfofile obj.type];
    cmd = [cmd ' -interp nearestneighbour'];
    
    obj.call_fsl(cmd);
end

%% Use Fslstats
%@param input image filepath+name
%@param fsl parameters
%@param outputfile to store result
function calcFslStats(obj, inputfile, param, outputfile)
    
    cmd = 'fslstats';
    cmd = [cmd ' ' inputfile];
    cmd = [cmd ' ' param '>' outputfile];
    
    obj.call_fsl(cmd);
end

%% Get center of gravity of an image
function cog = getCoG(obj, inputfile)
    
    cmd = 'fslstats';
    cmd = [cmd ' ' inputfile ' -c'];
    
    [status,output] = obj.call_fsl(cmd);
    cog = str2num(output);
end

%% Validdate T1 and AAL coregistation
% @param structural subjectpath
function fslViewAAL2T1(obj, subjectpath)
    %TODO parameterize in function
    transfo = 'AAL2T1';
    
    cmd = 'fslview';
    files = obj.getFilesIgnoreCase(subjectpath, obj.filepattern);
    aal = [subjectpath '/ROI_' transfo '.nii'];
    cmd = [cmd ' ' files ' '  aal];
    
    obj.call_fsl(cmd);
    
end

%% Validdate T1 and AAL coregistation
% @param structural subjectpath
% @param name of tranformation
function fslViewTransfo2T1(obj, subjectpath, transfo)
 
    cmd = 'fslview';
    files = obj.getFilesIgnoreCase(subjectpath, obj.filepattern);
    aal = [subjectpath '/ROI_' transfo '.nii'];
    cmd = [cmd ' ' files ' '  aal ' -l Yellow'];
    
    obj.call_fsl(cmd);
    
end


%% Validdate T2 and AAL coregistation
% @param structural subjectpath
function fslViewAAL2T2(obj, roipath, fMRIpath, fMRIpattern)
 
    cmd = 'fslview';
    files = spm_select('FPlist', fMRIpath, fMRIpattern);
    fMRI = cellstr(files);
    aal = [roipath '/r1ROI_AAL2T1.nii'];
    cmd = [cmd ' ' fMRI{1} ' '  aal];
    
    obj.call_fsl(cmd);
    
end

%% Validdate T2 and AAL coregistation
% @param structural subjectpath
function fslViewTransfo2T2(obj, roipath, fMRIpath, fMRIpattern, transfo)
 
    cmd = 'fslview';
    files = spm_select('FPlist', fMRIpath, fMRIpattern);
    fMRI = cellstr(files);
    aal = [roipath '/r1ROI_' transfo '.nii'];
    cmd = [cmd ' ' fMRI{1} ' '  aal];
    
    obj.call_fsl(cmd);
    
end

%% View overlayed images
% @param .img/.nii-files
function fslViewGeneric(obj, files)
 
    cmd = 'fslview';
    filecell = cellstr(files);

    for i = 1:length(filecell)
        cmd = [cmd ' ' filecell{i} ' '];
    end
    
    obj.call_fsl(cmd);
    
end

%% View overlayed images
% @param .img/.nii-files
% @param list of colors per file e.g.  Greyscale, Red-Yellow, Blue-Lightblue,
%               Red, Green, Blue, Yellow, Pink, Hot, Cool, Copper
function fslViewGenericColored(obj, files, colors)
 
    cmd = 'fslview';
    filecell = cellstr(files);

    for i = 1:length(filecell)
        cmd = [cmd ' ' filecell{i} ' -l ' colors{i}];
    end
    
    obj.call_fsl(cmd);
    
end

%% Calculates BET "robust" brain centre estimation
% @param mtpath
% @param mtpattern
% writes mtbet.nii into mt folder
function calcBET(obj, mtpath, mtpattern)
    
    cmd= 'bet';
    files = spm_select('FPlist', mtpath, mtpattern);
    mt = cellstr(files);
    betout = [mtpath '/mtbet.nii'];
    cmd = [cmd ' ' mt{1} ' ' betout ' -R'];
    
    obj.call_fsl(cmd);
    
end

%% Calculates BET "robust" brain centre estimation
% @param mtpath
% @param mtpattern
% writes mtbet.nii into mt folder
function calcBETGeneric(obj, path, inFile)
    
    cmd= 'bet';
    cmd = [cmd ' ' path inFile '.nii ' path inFile '_bet' obj.type ' -R'];
    
    obj.call_fsl(cmd);
    
end

%% split 4D images into 3D
% @param path
% @param pattern
% @param prefix
function split4D(obj, path, pattern, prefix)
    
    cd(path);
    cmd = 'fslsplit';
    files = spm_select('FPlist', path, pattern);
    files = cellstr(files);
    cmd = [cmd ' ' files{1} ' ' prefix];
    
    obj.call_fsl(cmd);
end

%% summarize images to 4d
% @param path
% @param pattern
% @param outputname
function make4d(obj, path, pattern, outputname)
    
    cd(path);
    cmd = 'fslmerge -a';
    files = spm_select('FPlist', path, pattern);
    files = cellstr(files);
    filestr = '';
    for i = 1:length(files)
        filestr = [filestr files{i} ' '];
    end
    cmd = [cmd ' ' outputname ' ' filestr];
    obj.call_fsl(cmd);
    
end


%% Calc eddy current correction
% @param path
% @param inFile
function eddy(obj, path, inFile)
    
%     cd(path);
    cmd = 'eddy_correct';
    
    in = [path inFile obj.type];
    out = [path inFile '_eddy' obj.type];
    
    cmd = [cmd ' ' in ' ' out ' 0']; %reference vol = 0
    
    obj.call_fsl(cmd);
    
end

%% Call fsl out of Matlab using the system command
% copied from fsl 5.0 /matlab folder
% Linux requires special commenting
% Currentversion tested for Ubuntu 13+12/Fsl5.0
% TODO eliminate hardcoded paths
function [status,output] = call_fsl(obj, cmd)
    % [status, output] = call_fsl(cmd)
    % 
    % Wrapper around calls to FSL binaries
    % clears LD_LIBRARY_PATH and ensures
    % the FSL envrionment variables have been
    % set up
    % Debian/Ubuntu users should uncomment as
    % indicated

%     fsldir=getenv('FSLDIR');
    fsldir = '/usr/share/fsl/5.0';

    % Debian/Ubuntu - uncomment the following
    fsllibdir=sprintf('%s/%s', fsldir, 'bin');

    if ismac
      dylibpath=getenv('DYLD_LIBRARY_PATH');
      setenv('DYLD_LIBRARY_PATH');
    else
      ldlibpath=getenv('LD_LIBRARY_PATH');
%       setenv('LD_LIBRARY_PATH');
      % Debian/Ubuntu - uncomment the following
      setenv('LD_LIBRARY_PATH',fsllibdir);
%       ldlibpath = '/usr/share/fsl/5.0/bin';
    end

    display(['Executing: ' cmd]);
    command = sprintf('/bin/sh -c ''. /usr/share/fsl/5.0/etc/fslconf/fsl.sh; %s''', cmd);
%     command = sprintf('%s''', cmd);
    [status,output] = system(command);
    if status == 2
        error(output);
    else
        display(['Execution finished: ' output]);
    end

    if ismac
      setenv('DYLD_LIBRARY_PATH', dylibpath);
    else
      setenv('LD_LIBRARY_PATH', ldlibpath);
    end
end

function files = getFilesIgnoreCase(obj, dir, filepattern)
    
    files = spm_select('FPlist', dir, filepattern);
    %workaround case insensitive filenames TODO do this fine
    if isempty(files)
        files = spm_select('FPlist', dir, lower(filepattern));
    end

end
    
end

end

