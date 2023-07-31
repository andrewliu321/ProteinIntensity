%% check boundary directions

%% instead of imsane the boundary axis, just load the raw
AP_i = record_AP{position};
DV_i = record_DV{position};

AP_i= im2bw(AP_i);
DV_i= im2bw(DV_i);
%% organize boundary axis

pause(1)

if exist('APi_index','var') == 0

    APi_index = cell(1);
    DVi_index = cell(1);

end

%load MIP if doesnt (after clearing everything)
% if exist('MaxIntensity','var') == 0
% end


% AP Surface Detected
% MaxIntensity_surf_stain = max(surf_stain, [], 3);
APi_index = CheckPixel(AP_i,MaxIntensity,APi_index,position); %if error, need to imshow(AP_i) and make sure its one continuous line without holes or extra pixels by itself
w1 = waitforbuttonpress;
% DV Surface Detected
DVi_index = CheckPixel(DV_i,MaxIntensity,DVi_index,position);
w2 = waitforbuttonpress;

%% should we flip DV? only run if its in wrong direction

if w2==1
    DVi_index{position} = flipud(DVi_index{position});
    disp('Flipped AP Boundary')
end
%%
if w1==1
    APi_index{position}= flipud(APi_index{position});
    disp('Flipped DV Boundary')
end
%% find center on the line and the index of the position

%for max project

% for surface of interest

SOI_Center = intersect(APi_index{position},DVi_index{position},'rows');

%if there's no intersecting point
if isempty(SOI_Center) == 1
   APi_temp = imdilate(AP_i, strel('line',2,90));
   DVi_temp = imdilate(DV_i, strel('line',2,90));
   
   %find pixel list of dilated lines
    clear xy_DVi
    clear xy_APi
    [xy_APi(:,1), xy_APi(:,2)] = find(APi_temp == 1);
    [xy_DVi(:,1), xy_DVi(:,2)] = find(DVi_temp == 1);
    
    DVi_cent = intersect(xy_APi,DVi_index{position},'rows');
    APi_cent = intersect(xy_DVi,APi_index{position},'rows');

    [~,APi_Center] = ismember(APi_cent,APi_index{position},'rows');
    [~,DVi_Center] = ismember(DVi_cent,DVi_index{position},'rows');
else

    [~,APi_Center] = ismember(SOI_Center,APi_index{position},'rows');
    [~,DVi_Center] = ismember(SOI_Center,DVi_index{position},'rows');

end

% % measure along the line 
% 1st and 2nd column = x,y index in image
% 3rd column = average intensity of 11 x 11 box
% 4th column = position index from neg to 0 to pos

% in case theres 2 center points, take the smaller one
APi_Center=min(APi_Center(APi_Center>0));
DVi_Center=min(DVi_Center(DVi_Center>0));

APi_Center_record{position}= APi_Center;
DVi_Center_record{position}= DVi_Center;

%%

function [x_index] = CheckPixel(matrix_ax,MIP_Stain,x_index,position)
%% input the axis matrix (AP,DV, AP_i, DV_i)

[x_index1,~] = bwboundaries(matrix_ax,8,'noholes');

%get endpoints of the line
matrix_ax_end = bwmorph(matrix_ax,'endpoints');
%get the pixel list of the two end points
[pts(:,1), pts(:,2)] = find(matrix_ax_end == 1);
%find where these poitns lie within the pixel list of the lines
[~,points] = ismember(pts,x_index1{1},'rows');
%remake index with only that list
x_index1{1} = x_index1{1}(min(points):max(points),:);

% x_index{1} = unique(x_index{1},'stable','rows');

size_x=size(x_index1{1},1);


overlay2 = cat(3, matrix_ax, MIP_Stain, 0*MIP_Stain);
figure(10)
hold off
imshow(overlay2)
figure(10)
hold on
title({'click mouse to not flip, enter to flip direction','anterior to post (left to right)','dorsal to ventral (top to bottom)'},'FontSize',25)

for s = 1:2:size_x/2
    
    plot(x_index1{1}(s,2),x_index1{1}(s,1),'*')
    pause(0.00005)

end

x_index{position}=x_index1{1};
end



