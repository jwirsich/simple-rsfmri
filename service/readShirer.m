% get atlas desacription from Shirer-Atrlas  [Shirer 2012, Cerebral Cortex]
% download https://findlab.stanford.edu/functional_ROIs.html (90 rois)
% 17-02-2017 Jonathan Wirsich / UIUC, Beckman
% 23-03-2017 Jonathan Wirsich / UIUC, Beckman read hippocampus
function shirer = readShirer(path)

    rsn = getRSN();
    
    for i = 1:length(rsn)
        
        idx=double.empty;
        
        subpath = [path rsn(i).name];
        
        %search for folders
        dirs = dir(subpath);
        dirIndex = find([dirs.isdir]);
        %cut directories /. and /..
        dirIndex = dirIndex(3:end);
        
        for j = 1:length(dirIndex)
            dirName = dirs(dirIndex(j)).name;
            idx(j) = str2double(dirName);
        end
        
        rsn(i).idx = idx;
    end
    
    shirer = rsn;
    
end

function shirer = getRSNOriginal()

    shirer(1).name = 'anterior_Salience';
    shirer(2).name = 'Auditory';
    shirer(3).name = 'Basal_Ganglia';
    shirer(4).name = 'dorsal_DMN';
    shirer(5).name = 'high_Visual';
    shirer(6).name = 'Language';
    shirer(7).name = 'LECN';
    shirer(8).name = 'post_Salience';
    shirer(9).name = 'Precuneus';
    shirer(10).name = 'prim_Visual';
    shirer(11).name = 'RECN';
    shirer(12).name = 'Sensorimotor';
    shirer(13).name = 'ventral_DMN';
    shirer(14).name = 'Visuospatial';

end

function shirer = getRSN()

    shirer(1).name = 'anterior_Salience';
    shirer(2).name = 'Auditory';
    shirer(3).name = 'Basal_Ganglia';
    shirer(4).name = 'dorsal_DMN';
    shirer(5).name = 'high_Visual';
    shirer(6).name = 'Language';
    shirer(7).name = 'LECN';
    shirer(8).name = 'post_Salience';
    shirer(9).name = 'Precuneus';
    shirer(10).name = 'prim_Visual';
    shirer(11).name = 'RECN';
    shirer(12).name = 'Sensorimotor';
    shirer(13).name = 'ventral_DMN';
    shirer(14).name = 'Visuospatial';
    shirer(15).name = 'Hippocampus';

end
