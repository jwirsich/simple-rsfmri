% srubbing imported from Marseille connectivity pipeline
% imports dVAR and calculates FD-power based on [Power2012, NeuroImage]
% 17-02-2017 Jonathan Wirsich / Connectlab

function [reject, accept, param] = scrub(dvars_file, sess_dir)
        
        %load motion parameters
        files = spm_select('List',sess_dir, '^rp.*');
        motion_file = files;
        tmp = load([sess_dir motion_file]);

        %framewise displacement

        translations = abs(diff(tmp(:, 4:6)));
        rotations = abs(diff(tmp(:, 1:3)));

        %estimate rotation displacment by simple spherical head model (diameter 50mm)
        fd_power = sum(translations, 2) + (50*3.141/180)* sum(rotations, 2);
        %FD_power = np.sum(translations, axis = 1) + (50*3.141/180)*np.sum(rotations, axis =1)

        %FD is zero for the first time point
        fd_power = [0; fd_power];

        if ~strcmp(dvars_file, '')
            %load svar
            tmp = load([sess_dir dvars_file]);
            dvars = zeros(length(tmp.intense),1);
            for i = 2:length(tmp.intense)
                dvars(i) = sqrt(mean((tmp.intense{i}-tmp.intense{i-1}).^2));
            end
            % %deltaBOLD
            dvars = dvars/100;
        else
            %skip dvars
            dvars=zeros(length(fd_power),1);
        end

        %define FD threshold: CPAC threshold = 0.2 / Power2012 -> 0.5
        %FD strict = 0.2 [Zalesky 2014, PNAS]
        %define DVARS threshold 0.5%deltaBOLD

        %loop mean timeseries
        reject = double.empty;
        for t = 1:length(dvars)
            if fd_power(t)>0.5 || dvars(t)>0.5
                %reject also frame before and after peak
                reject = [reject t-1 t t+1];
            end
        end
        %remove doubles
        reject = unique(reject);
        %put kept volume inidces in another list
        accept=1:length(dvars);
        accept = setdiff(accept,reject);
        
        display(['Scrubbed ' num2str(length(reject)) ' volumes'])
        %keep at least 150 datapoints!
        if length(reject)>200
            display(['Attention too many data points to be scrubbed (' num2str(length(reject)) ')'])
        end

        %printout excluded frames (to load into brainwaver) with ratio
        %write brainwaver segments for to cut out

        figure('name', 'Scrubbed')
        plot(1:length(dvars), fd_power, 1:length(dvars), dvars);
        
        param.dvar = dvars;
        param.fd = fd_power;
        
    end