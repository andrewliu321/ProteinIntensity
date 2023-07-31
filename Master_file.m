%%%% note: to measure expressions patterns by first drawing boundaries of a
%%%% set of images - this tutorial only includes 1 sample, followed by
%%%% ImSAnE, which will output a z-Stack, to be opened in FIJI and
%%%% peripodial erased, and lastly the signals measured. 


folder = 'C:\Users\Andrew\Desktop\ProteinIntensity\ImageData\';
% this is where you're image files exist

% input some parameters:
pixel_per_um = 4.158; 
z_steps = 0.3; %distance between each z slice
Num_Channels=4;
%Number of channels in the z-stack, example has 4 but depends on your
pre = 'Example '; 
%this could be genotype, is part of the name of the image before the #
post = ' L3';  
%lable type of larva/pupa, is part of the name of the image, after the #
save_name = ' for imsane.tif';
%save a new file ready for ImSAnE


%chose which file to load:  
for position = 1 
    %usually perform on a group of samples, so might be position = [1:n]  

        % run set up if first time
        if position == 1

        % run set up file
        addpath('C:\Users\Andrew\Desktop\ProteinIntensity\bfmatlab\')
        % Need to download bio-formats - included in github 
        addpath('imsane_AndrewCopy')
        % Need to download ImSAnE - included in github

        addpath('C:\Users\Andrew\Desktop\ProteinIntensity\Scripts\')

        cd 'C:\Users\Andrew\Desktop\ProteinIntensity\'
        set(0,'DefaultFigureWindowStyle','docked')
        end


    load_image
    %this will load 'Example 1 L3.tif'    
    %default bit depth is uint8, might need to change based on your image
    %right now it's loading 4 channels into separate variables:
    %Fat, PH3, Stain, Ds
   

    
    draw_pouch_axis
    %draw pouch and axis
    
    check_boundary_direction
    %make sure the direction are always Ant -> Post and Dorsal -> Ventral
end
% 
save('measured_variables_example',...
    'record_AP','record_DV',...
    'area_pouch','pouch_record',...
    'folder','position',...
    'DVi_index','APi_index')%,'-append')

%first time saving, don't include '-append' but 2nd time and so on, need to
%add to saved matlab file.


%% imsane
% before running ImSAnE, go to the imsane folder and run setup file

for position = 1 %again if doing a group of samples do [1:n]
dataDir = fullfile(folder);
onionOpts = struct('nLayers', 15, 'layerDistance', 1, 'sigma', 3,'makeMIP','MIP');
dataName = [pre,num2str(position),post,save_name];
save_name2 = ' imsaned.tif';
savedir = fullfile(folder, ['imsane\imsane_',pre,num2str(position),post,save_name2]);

run_imsane

after_imsane

implay(surf_im_edit(:,:,1,:))
% implay(ecad_im_edit)

end


save('measured_variables_example',...
    'size_x','size_y','size_z',...
    'Proper_Dist_AP','Proper_Dist_DV',...
    'DVi_index','APi_index','-append')%,...


%% measure the axis
%  run this after erasing imsane'd images to 
for position = 1 %again, for a group [1:n]
    
    after_erase
    % this specific example doesn't have any peripodial signal to erase

    % 9 zslice total to sum project
    start_z2 = 4 ;                
    end_z2   = 12;                
    % usually default from 4:12 is good.

    record_start2(position) =  start_z2;
    record_end2(position) = end_z2;

        
    SUM_Fat = (sum(Fat_im(:,:,start_z2:end_z2),3));
    SUM_Ds = (sum(Ds_im(:,:,start_z2:end_z2),3));

    measure_axis % measure 2 channels

end

save('measured_variables_example',...
    'Collect_DVi_index_sum_Ds','Collect_APi_index_sum_Ds',...
    'Collect_DVi_index_sum_Fat','Collect_APi_index_sum_Fat',...
    'record_end2','record_start2',...
    '-append')

%% Graphing and Output

graphing_and_output


save('measured_variables_example',...
    'Collect_DVi_index_sum_Ds','Collect_APi_index_sum_Ds',...
    'Collect_DVi_index_sum_Fat','Collect_APi_index_sum_Fat',...
    '-append')




