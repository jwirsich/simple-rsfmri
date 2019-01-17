%% write atlas coordinates from center of gravity
% Jonathan Wirsich 18/10/2017 / Connectlab
% depends on SPMFslTools in multimodal-private

shirerPath = '/media/jwirsich/DATAPART1/git/simple-rsfmri/atlas/shirer/';
shirer = readShirer(shirerPath);

fslTools = SPMfslTools('', '', 'nii');
count = 0;

for j = 1:length(shirer)    
     display(['Extracting RSN - ' shirer(j).name])
     for it = 1:length(shirer(j).idx)
            count = count+1;
            maskfile = fullfile([shirerPath filesep shirer(j).name filesep num2str(it, '%02d') ...
                filesep num2str(it, '%d') '.nii']);
            
            cog = fslTools.getCoG(maskfile);
            
            lines{count} = [num2str(count) ' ' shirer(j).name '_' num2str(shirer(j).idx(it)) ' ' num2str(cog)];
     end
end

for i = 1:length(lines)
    disp(lines{i})
end
