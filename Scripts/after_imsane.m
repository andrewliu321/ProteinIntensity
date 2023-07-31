%% get erased pouch imsane results

savedir2 = [savedir,'\fields\'];
savedir_7aup = [savedir2,'data_layer_m7\xy_index\xy\'];
savedir_6aup = [savedir2,'data_layer_m6\xy_index\xy\'];
savedir_5aup = [savedir2,'data_layer_m5\xy_index\xy\'];
savedir_4aup = [savedir2,'data_layer_m4\xy_index\xy\'];
savedir_3aup = [savedir2,'data_layer_m3\xy_index\xy\'];
savedir_2aup = [savedir2,'data_layer_m2\xy_index\xy\'];
savedir_1aup = [savedir2,'data_layer_m1\xy_index\xy\'];
savedir_0a = [savedir2,'data\xy_index\xy\'];
savedir_1a = [savedir2,'data_layer_p1\xy_index\xy\'];
savedir_2a = [savedir2,'data_layer_p2\xy_index\xy\'];
savedir_3a = [savedir2,'data_layer_p3\xy_index\xy\'];
savedir_4a = [savedir2,'data_layer_p4\xy_index\xy\'];
savedir_5a = [savedir2,'data_layer_p5\xy_index\xy\'];
savedir_6a = [savedir2,'data_layer_p6\xy_index\xy\'];
savedir_7a = [savedir2,'data_layer_p7\xy_index\xy\'];

clear surf_im_edit
surf_im_edit(:,:,1,1) = imread([savedir_7aup,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,2,1) = imread([savedir_6aup,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,3,1) = imread([savedir_5aup,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,4,1) = imread([savedir_4aup,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,5,1) = imread([savedir_3aup,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,6,1) = imread([savedir_2aup,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,7,1) = imread([savedir_1aup,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,8,1) = imread([savedir_0a,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,9,1) = imread([savedir_1a,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,10,1) = imread([savedir_2a,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,11,1) = imread([savedir_3a,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,12,1) = imread([savedir_4a,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,13,1) = imread([savedir_5a,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,14,1) = imread([savedir_6a,'cmp_1_1_T0000.tif']);
surf_im_edit(:,:,15,1) = imread([savedir_7a,'cmp_1_1_T0000.tif']);

surf_im_edit(:,:,1,2) = imread([savedir_7aup,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,2,2) = imread([savedir_6aup,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,3,2) = imread([savedir_5aup,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,4,2) = imread([savedir_4aup,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,5,2) = imread([savedir_3aup,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,6,2) = imread([savedir_2aup,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,7,2) = imread([savedir_1aup,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,8,2) = imread([savedir_0a,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,9,2) = imread([savedir_1a,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,10,2) = imread([savedir_2a,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,11,2) = imread([savedir_3a,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,12,2) = imread([savedir_4a,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,13,2) = imread([savedir_5a,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,14,2) = imread([savedir_6a,'cmp_2_1_T0000.tif']);
surf_im_edit(:,:,15,2) = imread([savedir_7a,'cmp_2_1_T0000.tif']);



%% load imsane output and save into tif file to erase
surf_im_edit = permute(surf_im_edit,[1,2,4,3]);

bfsave(surf_im_edit, [folder,pre,num2str(position),post,save_name2], 'dimensionOrder', 'XYCZT');%, 'BigTiff', true)





%% Calculate Actual Distance

%for AP use imsane function to calculate distance between each point

length_of_ax = size(APi_index{position},1);
Proper_Dist_temp_AP=zeros(length_of_ax,1);


parfor (i = 1:length_of_ax-1)

Proper_Dist_temp_AP(i+1)=xp.SOI.properLength(1,fliplr(APi_index{position}(i:i+1,:)),'xy');

end
% calculate cumulative sum
Proper_Dist_AP{position}=cumsum(Proper_Dist_temp_AP);

% find center point
AP_zero = APi_Center_record{position}; %find(Collect_APi_index_sum{position}(:,4)==0);
normalizing_factor_AP = Proper_Dist_AP{position}(AP_zero);
Proper_Dist_AP{position}(:,2)=Proper_Dist_AP{position}(:,1)-normalizing_factor_AP;

% %% for DV

length_of_ax = size(DVi_index{position},1);
Proper_Dist_temp_DV=zeros(length_of_ax,1);
parfor (i = 1:length_of_ax-1)

Proper_Dist_temp_DV(i+1)=xp.SOI.properLength(1,fliplr(DVi_index{position}(i:i+1,:)),'xy');

end

% calculate cumulative sum
Proper_Dist_DV{position}=cumsum(Proper_Dist_temp_DV);

% find center point
DV_zero = DVi_Center_record{position}; %find(Collect_DVi_index_sum{position}(:,4)==0);
normalizing_factor_DV = Proper_Dist_DV{position}(DV_zero);
Proper_Dist_DV{position}(:,2)=Proper_Dist_DV{position}(:,1)-normalizing_factor_DV;




%% functions









