%% load image files


%load file
im=bfopen([folder,pre,num2str(position),post,'.tif']);

zslices=size(im{1},1)/Num_Channels;
%set zeros
Ds=zeros(horzcat(size(im{1}{1}),zslices),'uint8'); %depends on what bit depth your image is
Fat=Ds;
Stain =Ds;
PH3 = Ds;

    for Z=1:zslices
        Fat(:,:,Z)=im{1}{4*Z-3}; % the first channel is Fat-GFP
        PH3(:,:,Z)=im{1}{4*Z-2}; % the second channel if anti-PH3
        Stain(:,:,Z)=im{1}{4*Z-1}; % the third channel is staining of Wg and En
        Ds(:,:,Z)=im{1}{4*Z};   %the fourth channel is anti-Ds
    end



%% show image

MIP_Fat_raw = max(Fat, [], 3);
MIP_Fat_raw=imadjust(uint8(MIP_Fat_raw));
%max projection of Fat-GFP

MIP_stain_raw = uint8(max(Stain, [], 3));
MIP_stain_raw=imadjust(uint8(MIP_stain_raw));
%max projection of staining of Wg and En

figure(4)
close(gcf)
figure(4)
ax1=subplot(1,2,1);
imshow(MIP_Fat_raw)

ax2=subplot(1,2,2);
imshow(MIP_stain_raw)







