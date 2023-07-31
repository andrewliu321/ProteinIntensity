%% after imsane and erase:
%load erased image
im=bfopen([folder,pre,num2str(position),post,save_name2]);


zslices=size(im{1},1)/2;
%set zeros
Ds_im=zeros(horzcat(size(im{1}{1}),zslices),'uint16');
Fat_im=Ds_im;

    for Z=1:zslices
        Ds_im(:,:,Z)=im{1}{2*Z};
        Fat_im(:,:,Z)=im{1}{2*Z-1};
    end
    
    implay(imadjustn(Ds_im))


