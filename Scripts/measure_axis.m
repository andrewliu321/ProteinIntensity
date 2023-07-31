%% measures axis in sum projected imsane,erased & max projected 

if exist('Collect_DVi_index_sum_Fat','var') == 0

    
    Collect_DVi_index_sum_Fat = cell(1);
    Collect_APi_index_sum_Fat = cell(1);
    Collect_DVi_index_sum_Ds = cell(1);
    Collect_APi_index_sum_Ds = cell(1);

end


f17=figure(17);
clf
imshow(imadjust(uint16(SUM_Fat)))
hold on
[Collect_APi_index_sum_Fat,~] = MeasureIntensity(APi_index,'y',SUM_Fat,APi_Center,position,Collect_APi_index_sum_Fat,25,f17,Ds2,0,'Avg');

f19=figure(19);
clf
imshow(imadjust(uint16(SUM_Ds)))
hold on
[Collect_APi_index_sum_Ds,~] = MeasureIntensity(APi_index,'y',SUM_Ds,APi_Center,position,Collect_APi_index_sum_Ds,25,f19,Ds2,0,'Avg');

% %% Do for DV

f18=figure(18);
clf
imshow(imadjust(uint16(SUM_Fat)))
hold on
[Collect_DVi_index_sum_Fat,~] = MeasureIntensity(DVi_index,'y',SUM_Fat,DVi_Center,position,Collect_DVi_index_sum_Fat,25,f18,Ds2,0,'Avg');

f20=figure(20);
clf
imshow(imadjust(uint16(SUM_Ds)))
hold on
[Collect_DVi_index_sum_Ds,~] = MeasureIntensity(DVi_index,'y',SUM_Ds,DVi_Center,position,Collect_DVi_index_sum_Ds,25,f20,Ds2,0,'Avg');
%%

if not(exist([folder,'measured_intensity_new']))
    mkdir([folder,'measured_intensity_new'])
end


figure(17)
title(['FatGFP DV Boundary, Pos =',num2str(position),' Window Size = 50'])
saveas(gcf,[folder,'measured_intensity_new\Fat_GFP_DVnum_',num2str(position),'.tif'],'tiff')
figure(18)
title(['FatGFP AP Boundary, Pos =',num2str(position),' Window Size = 50'])
saveas(gcf,[folder,'measured_intensity_new\Fat_GFP_APnum_',num2str(position),'.tif'],'tiff')
figure(19)
title(['anti Ds DV Boundary, Pos =',num2str(position),' Window Size = 50'])
saveas(gcf,[folder,'measured_intensity_new\anti_Ds_DVnum_',num2str(position),'.tif'],'tiff')
figure(20)
title(['anti Ds AP Boundary, Pos =',num2str(position),' Window Size = 50'])
saveas(gcf,[folder,'measured_intensity_new\anti_Ds_APnum_',num2str(position),'.tif'],'tiff')



%replace the old index with new:
Collect_APi_index_sum_Fat{position}(:,4)=Proper_Dist_AP{position}(:,2);
Collect_DVi_index_sum_Fat{position}(:,4)=Proper_Dist_DV{position}(:,2);
Collect_APi_index_sum_Ds{position}(:,4)=Proper_Dist_AP{position}(:,2);
Collect_DVi_index_sum_Ds{position}(:,4)=Proper_Dist_DV{position}(:,2);

%% plot measurements

if not(exist([folder,'measured_intensity_new']))
    mkdir([folder,'measured_intensity_new'])
end


f21=figure(21);
clf
hold on
scatter(Collect_APi_index_sum_Fat{position}(:,4)/pixel_per_um,Collect_APi_index_sum_Fat{position}(:,3),10,'g') 
scatter(Collect_APi_index_sum_Ds{position}(:,4)/pixel_per_um,Collect_APi_index_sum_Ds{position}(:,3),10,'r')
title(['Imsane Sum Projected along DV Boundary Pos =',num2str(position)])
xlabel('Length in micron')
ylabel('Fluor Intensity')
legend('Fat GFP','anti Ds')

saveas(gcf,[folder,'measured_intensity_new\num_',num2str(position),'_measured_DV.tif'],'tiff')

f22=figure(22);
clf
hold on

scatter(Collect_DVi_index_sum_Fat{position}(:,4)/pixel_per_um,Collect_DVi_index_sum_Fat{position}(:,3),10,'g')
scatter(Collect_DVi_index_sum_Ds{position}(:,4)/pixel_per_um,Collect_DVi_index_sum_Ds{position}(:,3),10,'r')
title(['Imsane Sum Projected Ds GFP along AP Boundary Pos =',num2str(position)])
xlabel('Length in micron')
ylabel('Fluor Intensity')
legend('Fat GFP','anti Ds')

saveas(gcf,[folder,'measured_intensity_new\num_',num2str(position),'_emeasured_AP.tif'],'tiff')





function [Collect_x_index,BoundaryInZ] = MeasureIntensity(x_index,color,ImageIntensityMatrix,xy_center,position,Collect_x_index,windowsize,figurenum,RawDs,DrawZ,AvgOrSum)
%% insert color = '*r' and x_index from CheckPixel output, ImageIntensityMatrix = DS image

%how big of a window to average? (0=1x1, 1= 3x3, 2=5x5, 3=7x7, etc)
% x_window_size=5;
% y_window_size=5;

xy_center=min(xy_center(xy_center>0));

Collect_x_index{position}=x_index{position};

%make 0's into NaN's
ImageIntensityMatrix=double(ImageIntensityMatrix);
ImageIntensityMatrix(ImageIntensityMatrix==0)=NaN;

%pad image to have NaN's around the space, pad amount by max of window size
% ImageIntensityMatrix=padarray(ImageIntensityMatrix,[windowsize windowsize],NaN);



%set temporary variable to hold index & mean intensity
temp1=nan(size(x_index{position},1),1);
temp2=nan(size(x_index{position},1),1);
BoundaryInZ=nan(size(RawDs,3),size(x_index{position},1));
size_z=size(RawDs,3);

parfor count1 = 1:size(x_index{position},1)

%     disp(count1)

% %old way with the x by y box averaging around a pixel
%     xx = x_index{position}(count1,2)+x_window_size;
%     yy = x_index{position}(count1,1)+y_window_size;
%     Collect_x_index{position}(count1,3) = nanmean(nanmean(ImageIntensityMatrix(yy-y_window_size:yy+y_window_size,xx-x_window_size:xx+x_window_size)));

    % for each pixel, define the center:
    center=[x_index{position}(count1,2) x_index{position}(count1,1)];
    
   
    if count1 < 4 
        %if 1st and 2nd pixel, then use itself as left bound
        left = center; 
    else
        %or else, left pixel = count1-2 (2 pixels to the left)
        left = [x_index{position}(count1-3,2) x_index{position}(count1-3,1)];
    end
    
    if count1 > size(x_index{position},1) - 3
        %if last or 2nd til last pixel, then use itself as right bound
        right = center;
    else
        %or else right pixel = count+2
        right = [x_index{position}(count1+3,2) x_index{position}(count1+3,1)];
    end

    
    %calculate slope - inverse negative slope = normal
    rise = (right(2)-left(2));
    run = (right(1)-left(1));
%     slope=rise/run;
    slope_inv = -run/rise;

    % if slope is horizontal, then just define the end points
    if rise==0
        solx=[center(1) center(1)];
        soly=[center(2)+windowsize center(2)-windowsize];
    else
        %if slope isn't horizontal, then solve for end points using slope and 
        %distance equations
        x = sym('x');
        y = sym('y');
        eqns = [(y-center(2))/(x-center(1))==slope_inv, windowsize^2 == (x-center(1))^2 + (y-center(2))^2];
        vars = [x y];
        [solx, soly] = solve(eqns,vars);
    end

    %turn solution into round numbers
    solx=round(double(solx));
    soly=round(double(soly));
    
    %define the points based on solutions...
    point1(count1,:) = [(solx(1)) (soly(1))];
    point2(count1,:) = [(solx(2)) (soly(2))];


    
    %calculate mean intensity    
    if AvgOrSum == 'Avg'
    temp1(count1)= nanmean(improfile(ImageIntensityMatrix,solx,soly));
    elseif AvgOrSum == 'Sum'
    temp1(count1)= nansum(improfile(ImageIntensityMatrix,solx,soly));
    elseif AvgOrSum == 'Max'
    temp1(count1)= max(improfile(ImageIntensityMatrix,solx,soly));
    end
    %record index based on center
    temp2(count1)=count1-xy_center;
    
    %hold on
    %plot(count,x_index{1}(count,3),color)
    
    %make a new z-stack that shows the contour of the AP/DV axis in Z...
    if DrawZ==1 %only if =1
        for z = 1:size_z

            temp3 = nanmean(improfile(RawDs(:,:,z),solx,soly));
            BoundaryInZ(z,count1)=temp3

        end
    end
end

Collect_x_index{position}(:,3)=temp1;
Collect_x_index{position}(:,4)=temp2;




figurenum;
hold on
plot(point1(:,1), point1(:,2),color)
plot(point2(:,1), point2(:,2),color)


%     imwrite(imadjust(uint16(SUM_ecad)),[folder,'measured_intensity_new\Pos_',num2str(position),'_ecad_SUM.tif'])
%     imwrite(MIP_ecad,[folder,'measured_intensity_new\Pos_',num2str(position),'_ecad_Max.tif'])

% plot(Collect_x_index{position}(1:1:end,4),Collect_x_index{position}(1:1:end,3))
end

%%
function [Anterior_index, Posterior_index, Ant_fit, Pos_fit] = BestFitPlot(Collect_x_index,position,color,order,figurenum)

j=0;
k=0;

if position==1
Anterior_index=cell(1);
Posterior_index=cell(1);
Ant_fit=cell(1);
Pos_fit=cell(1);
end

%separate each into Ant/Pos or Dorsal/Ventral
for x = 1:length(Collect_x_index{position})
    if Collect_x_index{position}(x,4)<0
        j=j+1;
        Anterior_index{position}(j,:)=Collect_x_index{position}(x,:);
    else
        k=k+1;
        Posterior_index{position}(k,:)=Collect_x_index{position}(x,:);
    end
end

%delete NaN's
Anterior_index{position}(any(isnan(Anterior_index{position}),2),:) = [];
Posterior_index{position}(any(isnan(Posterior_index{position}),2),:) = [];


%fit the anterior and posterior sizes separately
[Ant_fit{position},gof_a] = fit(Anterior_index{position}(:,4),Anterior_index{position}(:,3),'exp1');
[Pos_fit{position},gof_p] = fit(Posterior_index{position}(:,4),Posterior_index{position}(:,3),'exp1');

% plot(Ant_fit{position},Anterior_index{position}(:,4),Anterior_index{position}(:,3))
% hold on
% plot(Pos_fit{position},Posterior_index{position}(:,4),Posterior_index{position}(:,3))

% Plotting with colors...
figurenum;
hold on
scatter(Collect_x_index{position}(:,4),Collect_x_index{position}(:,3),10,color)
plot(Pos_fit{position},color,[0 max(Posterior_index{position}(:,4))],[0 0],'w')
plot(Ant_fit{position},color,[min(Anterior_index{position}(:,4)) 0],[0 0],'w')

%write the r^2 value of the fit in x label
str = ['Ant R^2 = ',num2str(round(gof_a.rsquare,2)),' ; Pos R^2 = ',num2str(round(gof_p.rsquare,2))];
if order == 1
    dim = [0.9 0.1 .5 .6];
elseif order ==2
    dim = [0.9 0.15 .5 .6];
elseif order ==3
    dim = [0.9 0.2 .5 .6];
elseif order ==4
    dim = [0.9 0.1000 0.3875 0.1706];
elseif order ==5
    dim = [0.9    0.1500    0.3875    0.1706];
elseif order ==6
    dim = [0.9    0.2000    0.3875    0.1706];
end

annotation(figurenum,'textbox',dim,'String',str,'Color',color,'FitBoxToText','on')
xlabel('Pixels from Center')
ylabel('Relative Intensity')

end

function [Boundary] = DrawBoundary(BoundaryInZ)
%% input DV axis

f16=figure(16)
clf
imshow(imadjust(uint8(BoundaryInZ)))

[x1,y1]=getpts(f16);
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

Boundary(:,1) = [ypixel1];
Boundary(:,2) = [xpixel1];

%check plot
f16
hold on
plot(x1,y1,'b')
plot(xpixel1,ypixel1,'c')%,'LineWidth',4)

end
