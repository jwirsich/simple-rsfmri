% get Nuissance regressors and correct fMRI timecourse
% import from Sepideh connectivity scripts
% 17-02-2017 Jonathan Wirsich / UIUC, Beckman
function f_regress_nuissance(sess_dir, atlas_dir, regress_globmean)    
        
    %% load motion parameters
    files = spm_select('List',sess_dir, '^rp.*');
    motion_file = files;
    rp = load([sess_dir motion_file]);
    % rp = zscore(rp); % converts to standard score along columns. doesn't change residuals of regstats (just changes amplitude of betas)!
    % detrending would change a lot though...
    
    %get timeseries
    tmp = load([sess_dir 'timeseries_raw.mat']);
    ts = tmp.timeseries;
    
    %load marsbar regressors
    wm = importdata([sess_dir 'wm.txt']);
    csf = importdata([sess_dir 'csf.txt']);
    
    if regress_globmean == 1
        gm = load([sess_dir filesep 'globalmean']);
    end
    
    %load threshold
    dim_ts = size(ts);
    regsout = zeros(dim_ts);
    
    %iterate over all regions
    for i = 1:dim_ts(1)
        %% detrend signal time course
        ROI_signal = detrend(ts(i,:));

        % checks
%         if find(isnan(ROI_signal)), error(['There are NaNs in region' int2str(i)]); end
        if find(isnan(ROI_signal)), continue; end
%         figure(1); imagesc([compssig rp]);
    %     pause(0.5); pause % wait for user to continue

        %% multiple regression
        if regress_globmean == 1
            tmp = regstats(ROI_signal', [csf wm rp gm.timeseries'], 'linear','r');
        else
            tmp = regstats(ROI_signal', [csf wm rp], 'linear','r');
        end
        regsout(i,:) = tmp.r;
    end

    %% save file
    save([sess_dir 'timeseries_regressed.mat'], 'regsout');

end