shirer = readShirer('/home/sepideh/Documents/git/simple-rsfmri/atlas/shirer/');
   
count = 1;
for j = 1:length(shirer)
%     display(['Extracting RSN - ' shirer(j).name])
    for i = 1:length(shirer(j).idx)
        display([num2str(count) ' ' shirer(j).name ' ' num2str(shirer(j).idx(i))])
        count = count + 1;
    end
end