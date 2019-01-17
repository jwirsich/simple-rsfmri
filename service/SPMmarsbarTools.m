% ectract summarized Nuissance regressors of T*2 from 5mm Rois
% imported from Marseille connectivity project
% 17-02-2017 Jonathan Wirsich / Connectlab
classdef SPMmarsbarTools
    
    properties
        filepattern
        MARS
    end
    
    methods
        
    function obj = SPMmarsbarTools(fp)
        obj.filepattern = fp;
        
        %TODO proper toolboxsetup
        % read any necessary defaults
%         if ~mars_struct('isthere', obj.MARS, 'OPTIONS')
          loadf = 1;
          obj.MARS.OPTIONS = [];
%         else
%           loadf = 0;
%         end
        [mbdefs sourcestr] = mars_options('Defaults');
        obj.MARS.OPTIONS = mars_options('fill', obj.MARS.OPTIONS, mbdefs);
%         mars_options('put');
        if loadf
          fprintf('Loaded MarsBaR defaults from %s\n',sourcestr);
        end
    end
    
    function obj = extractRois(obj, subjectpath)
        
        fprintf('Definition of CSF and white matter for: %s\n', subjectpath);
        conf.dir = subjectpath;
        conf.filepattern = obj.filepattern;
        defineROIcenter('init', conf);
        reply = input('Press Enter to continue', 's');
       
        %TODO evaluate: 
        %if extraction is slow its better to do the manual step before
        load([subjectpath filesep 'mb_csfwm_roicoords.mat'], 'coords');
        obj.defineSphereCenter(subjectpath, 'csf', coords{1});
        obj.defineSphereCenter(subjectpath, 'wm', coords{2});
                
        %extract ROI Data
        
        %CSF_roi.mat
        %spm no
        %modality fMRI
        %N sessions
        %functional images
        %scaling raw data 0
        
        %save to file
        %LCR_mdata
                
    end
    
    function defineSphereCenter(obj, subjectpath, tag, sphere_center)
        %make sphere
        roi.centre = sphere_center;
        roi.radius = 5;
        sphere_roi = maroi_sphere(roi);
        %sphere_roi = maroi_sphere(struct('centre', sphere_center, ...
        %                    'radius', 5));
        sphere_roi = label(sphere_roi, tag);
                
        saveroi(sphere_roi, fullfile(subjectpath, [tag '_roi.mat']));
        
%         if nargin < 2
          etype = 'default';
%         else
%           etype = varargin{2};
%         end
%         if nargin < 3
%           roi_list = '';
%         else
%           roi_list = varargin{3};
%         end

        varargout = {[]};

        [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Extract data', 0);

        % full options extraction

        marsD = [];
        %use no spm design
        [VY row] = SPMmars_image_scaling(marsD, subjectpath, obj.filepattern);

        % Summary function
        sumfunc = obj.MARS.OPTIONS.statistics.sumfunc;

        roi_list = spm_select('FPlist', subjectpath, [tag '_roi.mat']);
        % ROI names to objects
        o = maroi('load_cell', roi_list);

        % Do data extraction
        marsY = get_marsy(o{:}, VY, sumfunc, 'v');
        marsY = block_rows(marsY, row);
        if ~n_regions(marsY)
          msgbox('No data returned','Data extraction', 'warn');
          return
        end

        %set up armoire structure
%         o = marmoire;
%         spm_design_filter = mars_veropts('design_filter_spec');
%         filter_specs  = {[spm_design_filter(1,:);...
%         {'*_mdes.mat','MarsBaR: *_mdes.mat'}; ...
%             spm_design_filter(2:end,:)], ...
%         {'*_mdata.mat','MarsBaR data file (*_mdata.mat)'},...
%         {'*_mres.mat', 'MarsBaR results (*_mres.mat)'}};
%     
%         o = add_if_absent(o, 'def_design', ...
%             struct('default_file_name', 'untitled_mdes.mat',...      
%             'filter_spec', {filter_specs{1}},...
%             'title', 'Default design',...
%             'set_action','mars_arm_call(''set_design'',o,item,old_o)'));
%         o = add_if_absent(o, 'roi_data',...
%             struct('default_file_name', 'untitled_mdata.mat',...
%             'filter_spec', {filter_specs{2}},...
%             'title', 'ROI data',...
%             'set_action','mars_arm_call(''set_data'',o,item,old_o)'));
%         o = add_if_absent(o, 'est_design',...
%             struct('default_file_name', 'untitled_mres.mat',...
%             'filter_spec', {filter_specs{3}},...
%             'title', 'MarsBaR estimated design',...
%             'set_action', 'mars_arm_call(''set_results'',o,item,old_o)'));
    
        % set into armoire, and display
%         set_item_data(o, 'roi_data', marsY);
%         marsstruct = marsY.y_struct;
%         marsY.save_struct
            
          %get ROI struct
          marsstruct = y_struct(marsY);
        
%         mars_arm('set', 'roi_data', marsY);
%         mars_arm('show_summary', 'roi_data');

%         varargout = {marsY};
        
        %['mars_arm(''save_ui'', ''roi_data'', ' fw_st ');']       
%         fw_st = struct('force', 1, 'warn_empty', 1);
%         filename = struct([tag '_mdata'],  [tag '_mdata']);
        %['mars_arm(''save_ui'', ''est_design'', ' fw_st ');']
        %[varargout{1} o] = save_item_data_ui(o, varargin{:});
%         save_item_data(o, 'roi_data', fw_st, [tag '_mdata.mat']);
%         save_item_data(o, 'roi_data', fw_st);
        
        %TODO save .txt-file to disk
        
        fileID = fopen([subjectpath filesep tag '.txt'],'w');
        
        for i = 1:length(marsstruct.Y)
           fprintf(fileID,'%s\n', num2str(marsstruct.Y(i)) ); 
        end
        
        fclose(fileID);
            
    end
        
    end
    
end

