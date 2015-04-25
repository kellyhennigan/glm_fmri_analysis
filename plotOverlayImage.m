function [imgRgbs,olMasks,olVals,h,acpcSlices] = plotOverlayImage(nii,t1,cmap,c_range,plane,acpcSlices,doPlot)
%cmap,c_range,plane,acpcSlices)
% overlay rgb maps onto gray scaled anatomical background
%
% % INPUTS:
%
%     nii - nii file with vol = nii.data(:,:,:,1) used for overlay
%     t1 - anatomical nii file to be grayscaled background image
%     cmap - colormap for overlay; default is autumn
%     c_range - clip range for overlay;
%           abs(vol) < c_range(1) will be set to zero
%           abs(vol) > c_range(2) will be colored the same as c_range(2)
%     plane - 1 for sag, 2 for coronal, 3 for axial images; default is 3
%     acpcSlices - slices to plot in acpc coords; default is the mode
%           activation slice
%     doPlot - 0 to not plot figure; default is to plot (doPlot=1)
%
% % OUTPUTS:
%     imgRgbs - cell array of M x N x 3 arrays with r,g,b values along
%              the 3rd dimension
%     olMasks - " " binary arrays with 1s for voxels included in imgRgb
%              overlay and 0s otherwise
%     olVals - values plotted as an overlay. If overlay is in the same
%              space/has the dimesnions as the bg image, this is just the values
%              input in nii.data for the slice that is plotted.
%     h - cell array of figure handles
%     acpcSlices - acpcSlices that correspond to images & figures
%
%
% kjh March 2013
% based heavily on Bob's mrAnatOvelayMontage() function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol = nii.data(:,:,:,1); % overlay volume
xform = nii.qto_xyz;     % it's img to acpc xform
t1Vol = double(t1.data); % structural volume for background
t1Xform = t1.qto_xyz;    % it's img to acpc xform

    
% if an argument isn't given then assign default values
if notDefined('cmap')
    cmap = autumn;
end

if notDefined('c_range')
    c_range = [min(vol(vol~=0)), max(vol(vol~=0))];
end


if notDefined('plane')
    plane = 3; % default to axial slices
end

if notDefined('acpcSlices')
    [vi, vj, vk] = ind2sub(size(vol),find(vol));
    acpcCoords = mrAnatXformCoords(xform,[vi vj vk]);
    acpcSlices = round(unique(acpcCoords(:,plane)));
end

if size(acpcSlices,1)>1
    acpcSlices=acpcSlices'; % make acpcSlices a row vector
end

% plot overlay image? default is yes
if notDefined('doPlot')
    doPlot = 1;
end


lm=linspace(c_range(1),c_range(2),length(cmap))'; % values to color mapping

scSize = get(0,'ScreenSize'); % get screensize for pleasant figure viewing :)

%% threshold vol for overlay

% vol(abs(vol)<c_range(1))=0;
% vol(abs(vol)>c_range(2))=c_range(2);


%% autocrop borders

autoCropBorder = inf; % set to inf for no autocropping 
% 
if(isfinite(autoCropBorder))
    tmp = sum(t1Vol,3);
    x = find(sum(tmp,1)); x = [x(1) x(end)];
    y = find(sum(tmp,2)); y = [y(1) y(end)];
    tmp = squeeze(sum(t1Vol,1));
    z = find(sum(tmp,1)); z = [z(1) z(end)];
    clear tmp;
    pad = [-autoCropBorder autoCropBorder];
    x = x+pad; y = y+pad; z = z+pad;
    t1Vol = t1Vol(y(1):y(2),x(1):x(2),z(1):z(2));
    t1Xform = inv(inv(t1Xform)*[1 0 0 -y(1); 0 1 0 -x(1); 0 0 1 -z(1); 0 0 0 1]);
    clear tmp x y z pad
end



%% resample vol to be same size as t1Vol w/same dimensions

bb = t1Xform*[1,1,1,1;[size(t1Vol),1]]'; % bounding box for resampled vol
bb = bb(1:3,:)';
interp = [0 0 0 0 0 0]; % [1 1 1 0 0 0] for trilinear interpolation
% interp = [0 0 0 7 7 7]; % [1 1 1 0 0 0] for trilinear interpolation

vol = mrAnatResliceSpm(vol, inv(xform), bb, t1.pixdim, interp, false);

resNiiVol = vol; % to return as output
%% reorient t1Vol and vol depending on slice orientation

switch plane
    case 1 % sagittal; eyes looking left
        t1Vol = permute(t1Vol,[3 2 1]); t1Vol = flipdim(t1Vol,1);
        vol = permute(vol,[3 2 1]); vol = flipdim(vol,1);
        imgCoords = round(mrAnatXformCoords(inv(t1Xform),[acpcSlices;ones(1,length(acpcSlices));ones(1,length(acpcSlices))]));
        slStr = 'x = ';
        
    case 2 % coronal; top of the head is up
        t1Vol = permute(t1Vol,[3 1 2]); t1Vol = flipdim(t1Vol,1); 
        vol = permute(vol,[3 1 2]); vol = flipdim(vol,1); 
        imgCoords = round(mrAnatXformCoords(inv(t1Xform),[ones(1,length(acpcSlices));acpcSlices;ones(1,length(acpcSlices))]));
        slStr = 'y = ';
        
    case 3 % axial; eyes looking up
        t1Vol = permute(t1Vol,[2 1 3]);  t1Vol = flipdim(t1Vol,1);
        vol = permute(vol,[2 1 3]); vol = flipdim(vol,1);
        imgCoords = round(mrAnatXformCoords(inv(t1Xform),[ones(1,length(acpcSlices));ones(1,length(acpcSlices));acpcSlices]));
        slStr = 'z = ';
end
imgSlices = imgCoords(:,plane);


%% make rgb overlay images

for i = 1:length(imgSlices)
    
    % anatomical background image (grayscaled)
    t1Sl=t1Vol(:,:,imgSlices(i));
    t1Sl = t1Sl./max(max(t1Sl));  % scale vals from 0 to 1
    t1Rgb = repmat(t1Sl,[1,1,3]); % repping this means rgb values will be gray scaled
    
    % binary mask of overlay rgb image
    sl=vol(:,:,imgSlices(i));
    olMask = zeros(size(sl)); olMask(abs(sl)>0) = 1; % make a mask of for voxels w/overlay values
    olMask = repmat(olMask,[1,1,3]); % mask of rgb image
    idx_ol = find(olMask); % index for all vals included in rgb overlay
    
    % overlay rgb image w/rgb values from the colormap
    olRgb = olMask;
    [~,ci]=pdist2(lm,sl(abs(sl)>0),'euclidean','smallest',1); % ci is the colormap index for that voxel
    olRgb(idx_ol) = cmap(ci,:);   % fills in rgb values for each overlay voxel
    
    % rgb image w/anatomical background & rgb overlay
    imgRgb = t1Rgb;
    imgRgb(idx_ol) = olRgb(idx_ol);
    
    if doPlot
        h{i} = figure;
        pos = get(gcf,'Position');
        set(gcf,'Position',[scSize(3)-pos(3), scSize(4)-pos(4), pos(3), pos(4)]) % put the figure in the upper right corner of the screen
        image(imgRgb)
        axis equal; axis off;
        set(gca,'Position',[0,0,1,1]);
        text(size(sl,2)-20,size(sl,1)-20,[slStr,num2str(acpcSlices(i))],'color',[1 1 1])
    
     else
        h = []; % return empty figure handle if not plotting
   
    end
    
    
    % put output variables into a cell array 
    imgRgbs{i} = imgRgb; % rgb images for each slice plotted 
    olMasks{i} = olMask; % binary masks indicating voxels w/overlay values
    olVals{i} = sl;      % values corresponding to plotted overlays. If the overlay wasn't resampled to match the bg image, these will be the same as what was given as input nii.data(:,:,:,1)
    
    
end % slices

%% plot a colorbar as a separate figure
% 
% figure
% pos = get(gcf,'Position');
% set(gcf,'Position',[scSize(3)-pos(3), 60, pos(3), 100])
% image(1:length(lm))
% colormap(cmap)
% xticks = get(gca,'XTick');
% xticks = interp1([1 length(cmap)],c_range,xticks);
% xticks=num2str(xticks',2);
% set(gca,'XTickLabel',xticks)
% set(gca,'YTickLabel',[])
% set(gca,'box','off');
% set(gcf,'Color','w','InvertHardCopy','off','PaperPositionMode','auto');
% title('Colorbar')
% 
end % end of function 


