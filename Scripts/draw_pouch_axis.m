%% this file is to prepare images for ImSANE



%% height map of max project of Ds
% 
    [mip_Ds,h_map] = max(Fat, [], 3);


%% input wing pouch area
figure(4)
ax1=subplot(1,2,1);
title('Input Pouch Boundary')
[x2,y2]=getpts(ax1);
%round the points
x2=round(x2);
y2=round(y2);
%set zero to collect pixel list
xpixel2 = [];
ypixel2 = [];
    
for count = 1:size(x2,1)

    %if last point, connect to the first point
    if count == size(x2,1)
        count1=1;
    else %else, compare to the next point
        count1=count+1;
    end
    
    %x or y has bigger difference?
    dif_x2=abs(x2(count)-x2(count1));
    dif_y2=abs(y2(count)-y2(count1));

    if dif_x2 > dif_y2
        %which have bigger difference, find all points in between
        if x2(count)<x2(count1)
            xpixel2_temp=(x2(count):x2(count1));
        else
            xpixel2_temp=(x2(count):-1:x2(count1));
        end
        ypixel2_temp=round(interp1([x2(count);x2(count1)],[y2(count);y2(count1)],xpixel2_temp));
    else
        if y2(count)<y2(count1)
            ypixel2_temp=(y2(count):y2(count1));
        else
            ypixel2_temp=(y2(count):-1:y2(count1));
        end
        xpixel2_temp=round(interp1([y2(count);y2(count1)],[x2(count);x2(count1)],ypixel2_temp));
    end

    %add each segment to the previous one
    xpixel2 = [xpixel2 xpixel2_temp];
    ypixel2 = [ypixel2 ypixel2_temp]; 
end

%draw the circle that represents the pouch
pouch=zeros(size(h_map),'double');
indexes2 = sub2ind(size(h_map), ypixel2, xpixel2);
pouch(indexes2) = 1;
pouch = uint8(pouch);
pouch(pouch==1)=255;

%fill fill holes to get pouch
pouch1 = imfill(pouch,'holes');

pouch_record{position} = pouch1;
%records a mask for the pouch

area_pouch(position) = bwarea(pouch1);
%records area of the pouch

%% delete non pouch pixels in all the z slices
Stain2 = zeros(size(Ds),'uint8');
Ds2 = zeros(size(Ds),'uint16');
Fat2 = zeros(size(Ds),'uint8');
PH32 = zeros(size(Ds),'uint16');

num_of_z = size(Ds,3);
for zslice = 1:num_of_z
    temp = Stain(:,:,zslice);
    temp(pouch1==0)=NaN; 
    Stain2(:,:,zslice) = temp;
    
    temp = Ds(:,:,zslice);
    temp(pouch1==0)=NaN; 
    Ds2(:,:,zslice) = temp;
        
    temp = Fat(:,:,zslice);
    temp(pouch1==0)=NaN; 
    Fat2(:,:,zslice) = temp;
    temp = PH3(:,:,zslice);
    temp(pouch1==0)=NaN; 
    PH32(:,:,zslice) = temp;
    
end

%% make staining to draw axis on

% max project
% MaxIntensity = max((Stain2-Ds2*10), [], 3);
MaxIntensity = max(Stain2(:,:,:), [], 3);
MaxIntensity2 = max(Stain2(:,:,1:round(num_of_z)/3), [], 3); 
%decide to divide z by 3 or not
MIP_Fat = max(Fat2, [], 3);

% Adaptive Thresholding 
T=adaptthresh(MaxIntensity,0.6);
BW=imbinarize(MaxIntensity,T); 
BW1 = imerode(BW,strel('disk',1));
BW2 = imdilate(BW1,strel('disk',3));

overlay = cat(3, imadjust(MaxIntensity), BW2*125, imadjust(MaxIntensity2)*0);

%show partially projected stack to easier identification of WG
figure(4)
ax1=subplot(1,2,1);
imshow(MaxIntensity)
title('Draw Wg boundary (AP axis)')
% 
% if Num_Channels ==3;
% imshow(MIP_Ds)
% end
ax2=subplot(1,2,2);
imshow(overlay)

% %% make video to go through changes
% vid(:,:,1)=MaxIntensity;
% temp=zeros(1024,1024,'uint8');
% temp(BW==1)=255;
% vid(:,:,2)=temp;
% temp=zeros(1024,1024,'uint8');
% temp(BW1==1)=255;
% vid(:,:,3)=temp;
% 
% implay(vid)


%% input AP axis
[x,y]=getpts(ax1);
%round the points
x=round(x);
y=round(y);
%set zero to collect pixel list
xpixel = [];
ypixel = [];

for count = 1:size(x,1)-1

    %x or y has bigger difference?
    dif_x=abs(x(count)-x(count+1));
    dif_y=abs(y(count)-y(count+1));

    if dif_x > dif_y
        %which have bigger difference, find all points in between
        if x(count)<x(count+1)
            xpixel_temp=(x(count):x(count+1));
        else
            xpixel_temp=(x(count):-1:x(count+1));
        end
        ypixel_temp=round(interp1(x(count:count+1),y(count:count+1),xpixel_temp));
    else
        if y(count)<y(count+1)
            ypixel_temp=(y(count):y(count+1));
        else
            ypixel_temp=(y(count):-1:y(count+1));
        end
        xpixel_temp=round(interp1(y(count:count+1),x(count:count+1),ypixel_temp));
    end

    %add each segment to the previous one
    xpixel = [xpixel xpixel_temp];
    ypixel = [ypixel ypixel_temp]; 
end

%% input DV axis

%show whole project stack
figure(4)
ax1=subplot(1,2,1);
MaxIntensity3 = max(Stain2(:,:,round(num_of_z)*1/3:num_of_z), [], 3);
imshow(MaxIntensity3)
title('Draw En boundary (DV axis)')

[x1,y1]=getpts(ax1);
%round the points
x1=round(x1);
y1=round(y1);
%set zero to collect pixel list
xpixel1 = [];
ypixel1 = [];
    
for count = 1:size(x1,1)-1

    %x or y has bigger difference?
    dif_x1=abs(x1(count)-x1(count+1));
    dif_y1=abs(y1(count)-y1(count+1));

    if dif_x1 > dif_y1
        %which have bigger difference, find all points in between
        if x1(count)<x1(count+1)
            xpixel1_temp=(x1(count):x1(count+1));
        else
            xpixel1_temp=(x1(count):-1:x1(count+1));
        end
        ypixel1_temp=round(interp1(x1(count:count+1),y1(count:count+1),xpixel1_temp));
    else
        if y1(count)<y1(count+1)
            ypixel1_temp=(y1(count):y1(count+1));
        else
            ypixel1_temp=(y1(count):-1:y1(count+1));
        end
        xpixel1_temp=round(interp1(y1(count:count+1),x1(count:count+1),ypixel1_temp));
    end

    %add each segment to the previous one
    xpixel1 = [xpixel1 xpixel1_temp];
    ypixel1 = [ypixel1 ypixel1_temp]; 
end



%% check plot
ax1=subplot(1,2,1);
imshow(MaxIntensity)
hold on
plot(x,y,'r')
plot(xpixel,ypixel,'g')%,'LineWidth',4)
plot(x1,y1,'b')
plot(xpixel1,ypixel1,'c')%,'LineWidth',4)

ax1=subplot(1,2,2);
imshow(MIP_Fat)
hold on
plot(x,y,'r')
plot(xpixel,ypixel,'g')
plot(x1,y1,'b')
plot(xpixel1,ypixel1,'c')

%% make new matrix using the pixel list generated on AP and DV axis

% make image of AP line using the pixel list created earlier
AP=zeros(size(MaxIntensity),'uint8');
indexes = sub2ind(size(MaxIntensity), ypixel, xpixel);
AP(indexes) = 255;
% outside the boundary, erase the ap line
AP(pouch1==0)=0;

indexes1 = sub2ind(size(MaxIntensity), ypixel1, xpixel1);
DV = zeros(size(MaxIntensity),'uint8');
DV(indexes1)=255;
% outside the boundary, erase the ap line
DV(pouch1==0)=0;

%record AP and DV 
record_AP{position}=AP;
record_DV{position}=DV;



%% write new files

% % select new z range !!! DONE PREVIOUSLY IN FIRST SECTION
% start_z = 2;
% end_z = 14;

NewIm = zeros([size(Ds,1),size(Ds,2),10+size(Ds,3),2],'uint8'); 
%depends on bit image you have.


    NewIm(:,:,6:5+size(Ds,3),1) = Fat2(:,:,:);
    NewIm(:,:,6:5+size(Ds,3),2) = Ds2(:,:,:); 
 
    NewIm = permute(NewIm,[1,2,4,3]);

    if exist([folder,num2str(position),save_name])==2
      delete([folder,num2str(position),save_name]);
    end
    
    bfsave(NewIm, [folder,pre,num2str(position),post,save_name], 'dimensionOrder', 'XYCZT');

if not(exist([folder,'MIP_Stain_new']))
    mkdir([folder,'MIP_Stain_new'])
end

imwrite(MaxIntensity,[folder,'MIP_Stain_new\',num2str(position),'.tif'])



size_x(position) = size(NewIm,1);
size_y(position) = size(NewIm,2);
size_z(position) = size(NewIm,4);
























