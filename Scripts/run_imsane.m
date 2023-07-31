%% ImSAnE Wing Disc Tutorial
%
% In this tutorial we detect the apical surface of the disc proper cells
% and the peripodial cells on one side of a wing disc two hours after
% puparium formation. We read out actin and E-cadherin on these surfaces.
%
% This an example of how to detect and fit what we call planar surfaces,
% where the surface can be described as a height for every x,y.
% It also demonstrates using multiple channels and making multiple
% surfaces. Finally, because finding folDs as regions of high curvature is
% part of detecting the peripodial surface, this example shows how to
% calculate the surface metric and curvature.
%
% Note that ImSAnE is fully documented and additional information about
% available properties, methoDs and options can be found using the matlab
% documentation system. 
% 
% For example type:

%%

% doc surfaceDetection.planarDetector

%% Initialize the project
%
% We start by creating an experiment object, which holDs this metadata and 
% provides a frontend for a number of tasks such as loading the raw data.
% To construct the experiment object we need to pass dataDir, the directory 
% containing the raw data and projectDir, the directory where results of 
% the script will be saved.

% clear all; close all;

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

%dataDir = fullfile(scriptPath, 'rawData');
% dataDir = fullfile(folder);
projectDir = fullfile(scriptPath, 'projectFiles');

% test
xp = project.Experiment(projectDir, dataDir);

%%
% Next we set the metadata pertaining to the raw data files in the structure
% fileMeta. ImSAnE assumes that timepoints are saved individually and that 
% filenames in a timeseries are identical up to an integer specifying the 
% timepoint. Therefore we have
%
% * filenameFormat:               Filename, %u in the position of the integer
% * timePoints = [t1, t2, ..] :   List of times available. In this example 
% we have just a single time point 0.
% * stackResolution :             Stack resolution in micron.

um_per_pixel = 1/pixel_per_um;

fileMeta = struct();
fileMeta.dataDir = dataDir;
%fileMeta.filenameFormat = '032018_T%d.tif';
fileMeta.filenameFormat = dataName;
fileMeta.swapZT = 0;
fileMeta.stackResolution = [um_per_pixel um_per_pixel z_step]; %0.2405 0.2405 0.3
                                                 
fileMeta.rawStackSize = [size_x(position) size_y(position) 2 size_z(position)];
fileMeta.timePoints = [0];
%fileMeta.timePoints = [0:168];
fileMeta.nChannels = 2;

%% 
% In the structure expMeta we set general parameters for the surface
% analysis we will do. 
%
% * channelsUsed:   Which channels do we need.
% * channelColor:   Assign color to the channels, RGB = [1 2 3].
%                   In this example the first channel is E-cadherin, and 
%                   the second is actin. We want these in green and red,
%                   respectively.
% * dynamicSurface: Does the surface shape change with time?  
%                   For a single time point this is false. True is not yet
%                   supported.
% * jitterCorrection:   Not needed here.
% * detectorType:       Which type of surface detector will be used.
%                       We will look at only one side of the wing, which is
%                       a planar surface so we use planarDetector.
% * fitterType:         Which type of fitter will be used.
%                       We fit planar surfaces using Thin Plate Spline:
%                       tpsFitter.

expMeta = struct();
expMeta.description = 'wing disc';
expMeta.channelsUsed = [1 2];
expMeta.channelColor = [2 1];
expMeta.dynamicSurface = 0; %1
expMeta.jitterCorrection = 0;
% expMeta.fitTime = fileMeta.timePoints(1);
expMeta.detectorType = 'surfaceDetection.MIPDetector';
expMeta.fitterType = 'surfaceFitting.tpsFitter';

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

%%
% Finally we call initNew(), which reaDs the stack size from the first 
% available time point, then initializes fitter and detector and creates 
% fitOptions and detectOptions based on their defaults.

xp.initNew();

%% Load a time point from the data
%
% Now that we have set up the project, we can load a time point.
% loadTime sets xp.currentTime, loaDs the stack into xp.stack 
% and resets detector and fitter with the default options for that time.
%
% We rescale to unit aspect ratio to make detection work better later and
% visualize in the right proportions.

xp.loadTime(0);
% xp.rescaleStackToUnitAspect();

%% check the image content; 
figure(5)
 im = xp.stack.image.apply();
 imshow(max(im{1},[],3),[]);


% %%
% xp.stack is not an array but a Stack object.
% The easy way to look at a slice through the data is using getSlice.
% figure
% imshow(xp.stack.getSlice('z', 10), []);

%% Detect the surface
%
% planarDetector.detectSurface detects the surface as the position of the 
% maximal Gaussian z-derivative in some direction, i.e. the position of the
% largest intensity jump along some direction and smoothened over some
% scale.
%
% A number of detection options directly affect detection:
%
% * sigma :     Width of the Gaussian z-derivative.
% * channels :  Channels (summed) to use for detection.
% * zdir :      Dimension corresponding to z, minus flips direction.
% Flipping the direction can sometimes improve detection.
%
% Then there are options which filter the result and can be modified
% without redetecting:
%
% * maxIthresh:     Throw out points with MIP dimmer than this.
% * summedIthresh:  Throw out points with SIP dimmer than this.
% * sigZoutliers:   Remove height outliers after all other masks.
% * scaleZoutliers: Spatial scale of outlier removal.
%
% scaleZoutliers is the linear size of a region over which the
% distribution of height is computed, sigZoutliers is then a cutoff in
% units of standard deviation of this distribution to remove misdetected
% points far above or below the other points in the region.

%customOptions = struct('dilSize', 6, 'erSize', 20,'areaOpenSize', 1000, ...
%                        'fillSize', 1000);

detectOptions = struct('sigma', 5, 'channels', [1], 'zdir', 3,...
                        'maxIthresh', 0, 'summedIthresh', 0,...
                        'sigZoutliers', 10, 'scaleZoutliers', 100,...
                        'seedDistance', 20); %,'customOptions', customOptions); 
%'channels', [1 2 3]

% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.

xp.setDetectOptions(detectOptions);
xp.detectSurface();

%%% new section for other plots %%%
% %%
% % Different from the other detectors, the detected surface is
% % represented not only by a PointCloud object but also by an image
% % surfaceMatrix, containing z values for each xy.
% % Looking at this height map masked by the filters specified in
% % detectOptions one can judge how well the surface was detected.
% 
figure(6)
imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [],...
                                            'InitialMagnification', 100);

xp.detector.setManualMask(pouch_record{position})           
xp.detector.applyMasks
                                        % 
%                                         
% %%
% 
% imshow(mat2gray(xp.detector.surfaceMatrix))
% 
% bla = imextendedmax(xp.detector.surfaceMatrix, 10);
% imshow(bla, [])
% %%
% imshow(xp.detector.mask.*~bla.*xp.detector.surfaceMatrix, [])
% 
% %%
% imshow(cat(3, mat2gray(xp.detector.surfaceMatrix), imextendedmax(xp.detector.surfaceMatrix, 10),0*imregionalmax(xp.detector.surfaceMatrix)), [])
%%%%%%
% 
% %%
% % We can also inspect a point cloud cross section over the data with
% % detector.inspectQuality. In the pointCloud option, 'c' specifies the 
% % color cyan.
% figure
% subplot(2,1,1)
% inspectOptions= struct('dimension', 'x', 'value', 104, 'pointCloud', 'c');
% xp.detector.inspectQuality(inspectOptions, xp.stack);


%% 
% Or we can look at the point cloud in 3d, with some subsampling factor.
% ssfactor = 5;
% figure(4)
% xp.detector.pointCloud.inspect(ssfactor);

% %% Fit the surface for the disc proper cells KEVIN
% %
% % By detecting the largest intensity jump along z for each x,y in the
% % E-cad channel and filtering out local outliers we have found the apical
% % surface of the disc proper cells. We can now fit a smooth surface
% % representation to that.
% %
% % tpsFitter fits the pointcloud using a thin plate spline fit. It has the
% % following options:
% %
% % * griDsize:     Size of grid on which to generate fitted surface
% %               default [50 50], full size takes long.
% % * smoothing:    TPS smoothing parameter (default 1000).
% 
% fitmask = imfill(imdilate(xp.detector.mask, strel('disk',20)), 'holes');
% fitmask = imerode(fitmask,strel('disk',10));
% 
% fitOptions = struct('smoothing', 50, 'griDsize', [80 80], 'fitMask', fitmask);
% fitOptions.shift = 0;
% shift = xp.fitOptions().shift;
% 
% 
% xp.setFitOptions(fitOptions);
% xp.fitSurface();

%% Fit the surface for the disc proper cells TUTORIAL
%
% By detecting the largest intensity jump along z for each x,y in the
% E-cad channel and filtering out local outliers we have found the apical
% surface of the disc proper cells. We can now fit a smooth surface
% representation to that.
%
% tpsFitter fits the pointcloud using a thin plate spline fit. It has the
% following options:
%
% * griDsize:     Size of grid on which to generate fitted surface
%               default [50 50], full size takes long.
% * smoothing:    TPS smoothing parameter (default 1000).

fitOptions = struct('smoothing', 500, 'griDsize', [80 80]);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%
% We can visualize the result on a cross section with
% fitter.inspectQuality.

figure(3)
x=50;
subplot(2,1,1)
inspectOptions= struct('dimension', 'x', 'value', x, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

subplot(2,1,2)
inspectOptions= struct('dimension', 'x', 'value', x, 'pointCloud', 'c');%1150
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

%%
% % visualize 3D
% figure
% xp.fitter.inspectTPS()
% 

%%
% We now generate the Surface Of Interest. The charts to be generated are 
% specified in xp.fitter.charts. In this case there is only one, called
% 'xy'. 

xp.generateSOI();

%% Pull back the data to the surface
% 
% We pull back the data to the SOI using pullbackStack.

% xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);

% onionOpts = struct('nLayers', 5, 'layerDistance', 2, 'sigma', 3,'makeMIP','MIP');
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime, onionOpts);



%%
% To look at the pullback, we call the data field of the SOI at the right
% time, and get a particular patch from that with getPatch. A patch is a 
% part of a surface. In this case, there is only one called xy_index.
% Then we get the data in some patch in a particular coordinate system with
% getTransform. In this case there is only one coordinate system: xy.
% What we get is an object not only holding the image data but also
% metadata and methoDs to manipulate it. The actual data is obtained by
% calling the method apply. This returns a cell array with entries for each
% channel.

% xp.tIdx converts the time into an index in a list of time points
tidx = xp.tIdx(xp.currentTime);

% the first channel is Ds
channel = 1;

discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');

discProperImage = discProperPatch.getTransform('xy').apply{channel};

figure(7)
imshow(discProperImage, xp.stack.Ilim{1}, 'InitialMagnification', 100);

%% compare MIP and pullback


% 
% 
% 
% MIP = max(xp.stack.image.apply{1},[],3);
% blaMIP = mat2gray(MIP);
% blaSOI = mat2gray(discProperImage,double([min(MIP(:)) max(MIP(:))]));
% color = cat(3, blaSOI, blaMIP, 0*blaMIP);
% figure(8), imshow(color)
% 
% 
% % save to folder
% if not(exist([folder,image_type,'measured_intensity_new']))
%     mkdir([folder,image_type,'measured_intensity_new'])
% end
% 
% saveas(gcf,[folder,image_type,'measured_intensity_new\Pos_',num2str(position),'_image2.tif'],'tiff')

% imwrite(blaSOI, fullfile(dataDir, 'analysis', 'SOI.tif'));
% imwrite(blaMIP, fullfile(dataDir, 'analysis', 'MIP.tif'));
% imwrite(color, fullfile(dataDir, 'analysis', 'SOIMIP.tif'));
% 

%% Save the result
%
% Finally we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.

imwriteOptions = {'tif'};
% savedir = fullfile(folder, ['imsane\imsane_pos_',num2str(position)]);

options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         false);
xp.SOI.save(options)
