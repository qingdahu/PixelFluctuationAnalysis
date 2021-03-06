% March 24, 2020
% Please run the code section by seciton
% Make sure the required input file exists. mat files are generated by pixelFFT.m
% The filelocation is assumed to be in the same folder and file names match
% the .TIF files that were included with Hu et al. 2020.

%% Figure 1
clear 
close all
load('bead2-xosc.tif.mat')
figure(1)
imagesc(mean(flashbackstack,3))
axis equal ; colormap jet ; title([filename ' mean intensity'])
axis tight ; axis off ; colormap(gray)
pointsofinterest = [116,132;104,122;31,73;79,212;50,128]; %pointsofinterest = ginput;
timepoints = 0:0.064904:0.064904*399;
numofpts = size(pointsofinterest,1);
for points = 1:numofpts
    figure(1)
    hold on
    scatter( pointsofinterest(points,1),pointsofinterest(points,2),120,'x','red','LineWidth',4)
    hold off
    figure(2)
    subplot(numofpts,1,points)
    plot(timepoints,squeeze(flashbackstack(pointsofinterest(points,2),pointsofinterest(points,1),:)))
    axis tight ;
    ylim([0 3900]) ; set(gca,'FontSize',22)
    figure(3)
    subplot(5,1,points)
    plot(fq(2:200),squeeze(amplitudespectrums(pointsofinterest(points,2),pointsofinterest(points,1),2:200)))
    axis tight
    ylim([0 200]) ; set(gca,'FontSize',22)
    figure(2)
    subplot(5,1,5) ; xlabel('Time (s)')
    subplot(5,1,3) ; ylabel('Intensity (AU)')
    figure(3)
    subplot(5,1,5) ; xlabel('Frequency (1/s)')
    subplot(5,1,3) ; ylabel('Magnitude(AU)')
end

load('bead2-noosc.tif.mat')
timepoints = 0:0.064904:0.064904*399;
for points = 1:5
    figure(2)
    subplot(5,1,points)
    hold on
    plot(timepoints,squeeze(flashbackstack(pointsofinterest(points,2),pointsofinterest(points,1),:)))
    hold off
    figure(3)
    subplot(5,1,points)
    hold on
    plot(fq(2:200),squeeze(amplitudespectrums(pointsofinterest(points,2),pointsofinterest(points,1),2:200)))
    hold off
end

legend('x-oscillation','no oscillation')

%% figure 2

clear 
close all
load('bead2-xosc.tif.mat')
maxintensity=150;

figure(1) 
imagesc(amplitudespectrums(:,:,52))
axis tight; axis equal; axis off; colormap jet; %colorbar 
caxis([0 maxintensity]) ; %title([filename ' spectral amplitude at drive frequency'])

figure(2) 
imagesc(amplitudespectrums(:,:,51))
axis tight; axis equal; axis off; colormap jet; %colorbar 
caxis([0 maxintensity]) ; %title([filename ' spectral amplitude at fd-1'])

figure(3) 
imagesc(amplitudespectrums(:,:,53))
axis tight; axis equal; axis off; colormap jet; %colorbar 
caxis([0 maxintensity]) ; %title([filename ' spectral amplitude at fd+1'])

tempmean = mean(amplitudespectrums(:,:,[2:51,53:200]),3);
figure(4) 
imagesc(tempmean)
axis tight; axis equal; axis off; colormap jet; %colorbar 
caxis([0 maxintensity]) ; %title([filename ' mean'])

figure(5) 
imagesc(amplitudespectrums(:,:,52)-tempmean)
axis tight; axis equal; axis off; colormap jet; %colorbar 
caxis([0 maxintensity]) ; %title([filename ' meansubtract at fd'])

figure(6) 
imagesc(amplitudespectrums(:,:,51)-tempmean)
axis tight; axis equal; axis off; colormap jet; %colorbar 
caxis([0 maxintensity]) ; %title([filename ' meansubtract fd-'])

figure(7) 
imagesc(amplitudespectrums(:,:,53)-tempmean)
axis tight; axis equal; axis off; colormap jet; %colorbar 
caxis([0 maxintensity]) ;%title([filename ' mean substract fd+'])

intensitymask = squeeze(mean(flashbackstack,3)) ;
intensity_threshold = 1.2* mean(intensitymask(:)); % In practice this works well enough. A more precise method may be appropriate depending on the specific analysis.
intensitymask(intensitymask<intensity_threshold) = NaN;

tempamp = squeeze(amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3));
tempamp(isnan(intensitymask)) = NaN;
figure(8)
fig1=pcolor(tempamp);
set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse');rectangle('Position',[1 1 255 255])
axis tight ; axis equal  ; axis off
caxis([0 150]) ; colormap jet ; set(gca,'FontSize',22); %colorbar 
%title([filename ' original masked'])

tempamp = squeeze(amplitudespectrums(:,:,51) - mean(amplitudespectrums(:,:,[2:51,53:200]),3));
tempamp(isnan(intensitymask)) = NaN;
figure(9)
fig1=pcolor(tempamp);
set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse');rectangle('Position',[1 1 255 255])
axis tight ; axis equal  ; axis off
caxis([0 150]) ; colormap jet ; set(gca,'FontSize',22); %colorbar 
%title([filename ' original masked'])

load('/Volumes/PhD/botvinick/data/20191219/bead2-noosc.tif.mat')
maxintensity=150;
figure(10) 
imagesc(amplitudespectrums(:,:,52))
axis tight ; axis equal  ; axis off 
colormap jet ; caxis([0 maxintensity]) ; %colorbar 
%title([filename ' spectral amplitude at drive frequency'])

intensitymask = squeeze(mean(flashbackstack,3)) ;
intensity_threshold = 1.2* mean(intensitymask(:)); 
intensitymask(intensitymask<intensity_threshold) = NaN;

tempamp = squeeze(amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3));
tempamp(isnan(intensitymask)) = NaN;
figure(11)
fig1=pcolor(tempamp);
set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse');rectangle('Position',[1 1 255 255])
axis tight ; axis equal  ; axis off
caxis([0 150]) ; colormap jet ; set(gca,'FontSize',22); %colorbar 
%title([filename ' original masked'])



%% figure 3 

clear 
close all
load('bead2-noosc.tif.mat')
tempangles = angle(myffts(:,:,52));
figure(1) 
imagesc(tempangles)
axis tight ; axis equal  ; axis off  ; %colorbar 
colormap hsv
title([filename ' mean subtracted'])

figure(2)
histogram(tempangles(:),100)
set(gca,'FontSize',22)
xlim([-pi pi])
ylim([0 1000])

figure(3)
fullpha = angle(myffts(:,:,52));
fullmag = squeeze(amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3));
scatter(fullpha(:),fullmag(:),3,'r')
set(gca,'FontSize',22)
xlim([-pi pi])
ylim([0 400])


load('bead2-xosc.tif.mat')
tempangles = angle(myffts(:,:,52));
figure(4) 
imagesc(tempangles)
axis tight ; axis equal  ;axis off 
colormap hsv; %colorbar 

figure(4)
histogram(tempangles(:),100)
set(gca,'FontSize',22)
xlim([-pi pi])
ylim([0 1000])

figure(6)
fullpha = angle(myffts(:,:,52));
fullmag = squeeze(amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3));
scatter(fullpha(:),fullmag(:),3,'r')
set(gca,'FontSize',22)
xlim([-pi pi])
ylim([0 400])


intensitymask = squeeze(mean(flashbackstack,3)) ;
intensity_threshold = 1.2* mean(intensitymask(:));
intensitymask(intensitymask<intensity_threshold) = NaN;
tempangles = angle(myffts(:,:,52));
tempangles(isnan(intensitymask)) = NaN;

figure(7)
fig1=pcolor(tempangles);
set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse')  ;rectangle('Position',[1 1 255 255])
axis tight ;axis equal  ;axis off
colormap hsv ;set(gca,'FontSize',22); %colorbar 
%title([filename ' original masked'])

sze = 3;
[Ioop, Ineigh] = hoodOOP(squeeze(angle(myffts(:,:,52))),sze);
figure(8) 
IOOPMASKED = Ioop;
IOOPMASKED(isnan(intensitymask)) = nan;
fig1=pcolor(IOOPMASKED);
set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse')  ;rectangle('Position',[1 1 255 255])
axis tight ;axis equal  ;axis off
caxis([0 1]);colormap jet; colorbar ; set(gca,'FontSize',22)
%title(['masked neighborhood OOP with sze' num2str(sze)])

%%
% Ioop2 = Ioop(:);
% fullmag2= imfilter( fullmag,ones(7)/49);
% fullmag2= fullmag2(:);
% Ioop2(isnan(intensitymask(:))) = [];
% fullmag2(isnan(intensitymask(:))) = [];
% 
% 
% fiteq = 'a/(x-1)+a';
% [bestfit,goodness,output] = fit(Ioop2,fullmag2,fiteq,'StartPoint',[-1] );
% 
% figure(10)
% plot(bestfit,Ioop2,fullmag2)
% set(gca,'FontSize',28)
% xlim([0 1])
% ylim([0 150])
% legend off
% xlabel('') 
% ylabel('') 

%% figure 4

clear
close all

beadname = 'bead5_';
listfilenames = ["0mv.tif.mat" , ...
            "2mv.tif.mat" , ...
            "5mv.tif.mat"  , ...
            "7mv.tif.mat", ...
            "10mv.tif.mat", ...
            "12mv.tif.mat"];
        
for figurenum = 1:6
    load([beadname char(listfilenames(figurenum))])  
    
    intensitymask = squeeze(mean(flashbackstack,3)) ;
    intensity_threshold = 1.2* mean(intensitymask(:));
    intensitymask(intensitymask<intensity_threshold) = NaN;
    tempamp = squeeze(amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3));
    %tempamp(isnan(intensitymask)) = NaN;
    tempangles = angle(myffts(:,:,52));
    tempangles(isnan(intensitymask)) = NaN;
    
    figure(figurenum*3-2)
    fig1=pcolor(tempamp);
    set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse');rectangle('Position',[1 1 255 255])
    axis tight ; axis equal  ; axis off
    caxis([0 75]) ; colormap jet ; set(gca,'FontSize',22); %colorbar 

    figure(figurenum*3-1)
    fig1=pcolor(tempangles);
    set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse')  ;rectangle('Position',[1 1 255 255])
    axis tight ;axis equal  ;axis off
    colormap hsv ;set(gca,'FontSize',22); %colorbar 

   
%     figure(figurenum*3)
%     intensitymask = squeeze(mean(flashbackstack,3)) ;
%     intensity_threshold = 1.2* mean(intensitymask(:)); % In practice this works well enough. A more precise method may be appropriate depending on the specific analysis.
%     intensitymask(intensitymask<intensity_threshold) = NaN;
%     sze = 3;
%     [Ioop, Ineigh] = hoodOOP(squeeze(angle(myffts(:,:,52))),sze);
%     IOOPMASKED = Ioop;
%     IOOPMASKED(isnan(intensitymask)) = nan;
%     fig1=pcolor(IOOPMASKED); rectangle('Position',[1 1 255 255])
%     set(fig1, 'EdgeColor', 'none'); axis tight ;axis equal ;axis off ;caxis([0 1])
%     colormap jet; set(gca,'YDir','reverse') ; %set(gca,'FontSize',22); colorbar %title(['masked neighborhood OOP with sze' num2str(sze)])
%  
end


%% figure 5

clear
close all
listfilenames = [   "bead2-xosc.tif.mat" , ...
                    "bead2-yosc.tif.mat" , ...
                    "bead6_10mvxosc.tif.mat"  , ...
                    "bead6_10mvyosc.tif.mat", ...
                    "bead7_10mvxosc.tif.mat"  , ...
                    "bead7_10mvyosc.tif.mat", ...
                    "bead12_10mvxosc.tif.mat"  , ...
                    "bead12_10mvyosc.tif.mat"];
        

for figurenum = 1:8
    load( listfilenames(figurenum) )  

    intensitymask = squeeze(mean(flashbackstack,3)) ;
    intensity_threshold = 1.2* mean(intensitymask(:));
    intensitymask(intensitymask<intensity_threshold) = NaN;
    tempamp = squeeze(amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3));
    tempamp(isnan(intensitymask)) = NaN;
    tempangles = angle(myffts(:,:,52));
    tempangles(isnan(intensitymask)) = NaN;
    
    figure(figurenum*3-2)
    fig1=pcolor(tempamp);
    set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse');rectangle('Position',[1 1 255 255])
    axis tight ; axis equal  ; axis off
    caxis([0 75]) ; colormap jet ; set(gca,'FontSize',22); %colorbar 

    figure(figurenum*3-1)
    fig1=pcolor(tempangles);
    set(fig1, 'EdgeColor', 'none'); set(gca,'YDir','reverse')  ;rectangle('Position',[1 1 255 255])
    axis tight ;axis equal  ;axis off
    colormap hsv ;set(gca,'FontSize',22); %colorbar 

    figure(figurenum*3)
    intensitymask = squeeze(mean(flashbackstack,3)) ;
    intensity_threshold = 1.2* mean(intensitymask(:)); % In practice this works well enough. A more precise method may be appropriate depending on the specific analysis.
    intensitymask(intensitymask<intensity_threshold) = NaN;
    sze = 3;
    [Ioop, Ineigh] = hoodOOP(squeeze(angle(myffts(:,:,52))),sze);
    IOOPMASKED = Ioop;
    IOOPMASKED(isnan(intensitymask)) = nan;
    %fig1=pcolor(Ioop);
    fig1=pcolor(IOOPMASKED); rectangle('Position',[1 1 255 255])
    set(fig1, 'EdgeColor', 'none'); axis tight ;axis equal ;axis off ;caxis([0 1])
    colormap jet; set(gca,'YDir','reverse') ; %set(gca,'FontSize',22); colorbar 
    %title(['masked neighborhood OOP with sze' num2str(sze)])
   
end


%% figure 6
clear
close all
listfilenames = [   "bead2-xosc.tif.mat" , ...
                    "bead2-yosc.tif.mat" , ...
                    "bead6_10mvxosc.tif.mat"  , ...
                    "bead6_10mvyosc.tif.mat", ...
                    "bead7_10mvxosc.tif.mat"  , ...
                    "bead7_10mvyosc.tif.mat", ...
                    "bead12_10mvxosc.tif.mat"  , ...
                    "bead12_10mvyosc.tif.mat"];    
centerishs = [127 126; 127 129 ; 128 126 ;132 128];

intensitymaskstream = [];
xtempampstream      = [];
ytempampstream      = [];
anglestream         = [];
distancestream      = []; 

for figurenum = 1:4
    
    load( listfilenames(figurenum*2-1) )  %xosc
    intensitymask = squeeze(mean(flashbackstack,3)) ;
    intensity_threshold = 1.2* mean(intensitymask(:)); % In practice this works well enough. A more precise method may be appropriate depending on the specific analysis.
    intensitymask(intensitymask<intensity_threshold) = 0;
    intensitymask(intensitymask>intensity_threshold) = 1;
    
    xtempamp = amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3);
    xtempamp2 = xtempamp;
    xtempamp2(xtempamp<0) = 0;
    
    load( listfilenames(figurenum*2) )  %yosc
    ytempamp = amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3);
    ytempamp2 = ytempamp;
    ytempamp2(ytempamp<0) = 0;

    centerish = centerishs(figurenum,:);
    xmatrix = repmat(1:256, [256 1]);
    ymatrix = repmat([1:256]', [1 256]);
    distancematrix = sqrt( (xmatrix - centerish(1) ).^2 + (ymatrix - centerish(2) ).^2)*0.082;
    anglematrix = atan((xmatrix - centerish(1) )./(ymatrix - centerish(2) ));
    anglematrix(centerish(2),centerish(1)) = 0;
    
    %intensitymaskstream_ = 
    xtempampstream_                         =xtempamp2;
    xtempampstream_(intensitymask==0)       =NaN;
    xtempampstream_                         =xtempampstream_(:);
    xtempampstream_(isnan(xtempampstream_)) =[];
    
    ytempampstream_                        =ytempamp2;
    ytempampstream_(intensitymask==0)      =NaN;
    ytempampstream_                        =ytempampstream_(:);
    ytempampstream_(isnan(ytempampstream_))=[];
    
    anglestream_                            =anglematrix;
    anglestream_(intensitymask==0)          =NaN;
    anglestream_                            =anglestream_(:);
    anglestream_(isnan(anglestream_))       =[];
    
    distancestream_                         =distancematrix;
    distancestream_(intensitymask==0)       =NaN;
    distancestream_                         =distancestream_(:);
    distancestream_(isnan(distancestream_)) =[];
    
    %intensitymaskstream = [intensitymaskstream, intensitymaskstream_];
    xtempampstream      = [xtempampstream; xtempampstream_];
    ytempampstream     = [ytempampstream; ytempampstream_];
    anglestream         = [anglestream; anglestream_];
    distancestream      = [distancestream; distancestream_]; 
    
end

%anglestream_x = anglestream+pi/2;
%anglestream_x(anglestream_x>pi/2) = pi- anglestream_x(anglestream_x>pi/2);
%anglestream_y = abs(anglestream);

binnedmatrix_on = [];
binnedmatrix_off = [];
binnedmatrix_on2 = [];
binnedmatrix_off2 = [];
%binnedmatrix = [];
rawmatrix_on = {};
rawmatrix_off = {};

for j = 1:2:9
        binnedmatrix_off(j) = mean([xtempampstream( abs(anglestream)< pi/8 & abs(anglestream)>7*pi/8 & distancestream> j & distancestream< j+2 ) ; ytempampstream(abs(anglestream)>3*pi/8 & abs(anglestream)<5*pi/8 & distancestream> j & distancestream< j+2 )]); 
        binnedmatrix_on(j) = mean([ytempampstream(abs(anglestream)< pi/8 & abs(anglestream)>7*pi/8 & distancestream> j & distancestream< j+2 ); xtempampstream(abs(anglestream)>3*pi/8 & abs(anglestream)<5*pi/8 & distancestream> j & distancestream< j+2 )]);      
        binnedmatrix_off2(j) = median([xtempampstream( (abs(anglestream)< pi/8 | abs(anglestream)>7*pi/8) & distancestream> j & distancestream< j+2 ) ; ytempampstream(abs(anglestream)>3*pi/8 & abs(anglestream)<5*pi/8 & distancestream> j & distancestream< j+2 )]); 
        binnedmatrix_on2(j) = median([ytempampstream((abs(anglestream)< pi/8 | abs(anglestream)>7*pi/8) & distancestream> j & distancestream< j+2 ) ; xtempampstream(abs(anglestream)>3*pi/8 & abs(anglestream)<5*pi/8 & distancestream> j & distancestream< j+2 )]);      
        %binnedmatrix(i,j) = mean([xtempampstream(anglestream_x> i*pi/8-pi/8 & anglestream_x< i*pi/8 & distancestream> j & distancestream< j+2 ) ; %ytempampstream(anglestream_y> i*pi/8-pi/8 & anglestream_y< i*pi/8 & distancestream> j & distancestream< j+2) ]); 
        rawmatrix_off{j} = [xtempampstream( (abs(anglestream)< pi/8 | abs(anglestream)>7*pi/8) & distancestream> j & distancestream< j+2 ) ; ytempampstream(abs(anglestream)>3*pi/8 & abs(anglestream)<5*pi/8 & distancestream> j & distancestream< j+2 )]; 
        rawmatrix_on{j} = [ytempampstream((abs(anglestream)< pi/8 | abs(anglestream)>7*pi/8) & distancestream> j & distancestream< j+2 ); xtempampstream(abs(anglestream)>3*pi/8 & abs(anglestream)<5*pi/8 & distancestream> j & distancestream< j+2 )];  
end

figure(2)
plot(2:2:10,binnedmatrix_on2(1:2:9),'c',2:2:10,binnedmatrix_off2(1:2:9),'m','LineWidth',4);xlim([2 10]); ylim([0 9])
%plot(2:2:10,binnedmatrix_on2(1:2:9),'co',2:2:10,binnedmatrix_off2(1:2:9),'mo','LineWidth',4);xlim([2 10]); ylim([0 9])

% hold on
%plot(2:2:10,binnedmatrix_on(1:2:9),'c--',2:2:10,binnedmatrix_off(1:2:9),'m--','LineWidth',4); 

% scatter(1.8+0.2*rand(size(rawmatrix_on{1})), rawmatrix_on{1},1,'c'); xlim([1 11]); ylim([0 200])
% scatter(3.8+0.2*rand(size(rawmatrix_on{3})), rawmatrix_on{3},1,'c')
% scatter(5.8+0.2*rand(size(rawmatrix_on{5})), rawmatrix_on{5},1,'c')
% scatter(7.8+0.2*rand(size(rawmatrix_on{7})), rawmatrix_on{7},1,'c')
% scatter(9.8+0.2*rand(size(rawmatrix_on{9})), rawmatrix_on{9},1,'c')
% scatter(2+0.2*rand(size(rawmatrix_off{1})), rawmatrix_off{1},1,'m')
% scatter(4+0.2*rand(size(rawmatrix_off{3})), rawmatrix_off{3},1,'m')
% scatter(6+0.2*rand(size(rawmatrix_off{5})), rawmatrix_off{5},1,'m')
% scatter(8+0.2*rand(size(rawmatrix_off{7})), rawmatrix_off{7},1,'m')
% scatter(10+0.2*rand(size(rawmatrix_off{9})), rawmatrix_off{9},1,'m')

% hold off
%  0,1,1;1,1,1;1,0,1 xlabel('Distance away from center of bead (\mum)'); ylabel('')
set(gca,'FontSize',22)
%legend('Median parallel to oscillation','Median perpendicular to oscillation','Mean parallel to oscillation','Mean perpendicular to oscillation')
legend('Parallel to oscillation','Perpendicular to oscillation')

% statistical test for difference in median
for i = 1:2:9
   ranksum(rawmatrix_off{i},rawmatrix_on{i}) 
end
% for i = 1:2:9
%    ranksum(nonzeros(rawmatrix_off{i}),nonzeros(rawmatrix_on{i}) )
% end

figure(3) 
angleoverlay = zeros(256);
angleoverlay((abs(anglematrix)<pi/8 | abs(anglematrix)>7*pi/8)& distancematrix>1 & distancematrix<11 ) = 1;
angleoverlay( abs(anglematrix)>3*pi/8 & abs(anglematrix)<5*pi/8 & distancematrix>1 & distancematrix<11 ) = -1;
imagesc(angleoverlay); axis equal;axis tight;axis off
colormap([0,1,1;1,1,1;1,0,1]);caxis([-1 1]); %colorbar; 

figure(4)
imagesc(angleoverlay); axis equal;axis tight;axis off
colormap([1,0,1;1,1,1;0,1,1]);caxis([-1 1]); %colorbar; 

%%
                
for figurenum = 1:4
    
    load( listfilenames(figurenum*2-1) )  %xosc
    intensitymask = squeeze(mean(flashbackstack,3)) ;
    
%     figure(20+figurenum)
%     imagesc(intensitymask)
    
    intensity_threshold = 1.2* mean(intensitymask(:)); % In practice this works well enough. A more precise method may be appropriate depending on the specific analysis.
    intensitymask(intensitymask<intensity_threshold) = 0;
    intensitymask(intensitymask>intensity_threshold) = 1;
    xtempamp = amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3);
    xtempamp2 = xtempamp;
    xtempamp2(xtempamp<0) = 0;
    xtemptotamp = imgaussfilt( amplitudespectrums(:,:,[2:51,53:200])- mean(amplitudespectrums(:,:,[2:51,53:200]),3),2);
    xampthreshold = prctile(xtemptotamp(:),99.99);
    
    load( listfilenames(figurenum*2) )  %yosc
    ytempamp = amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3);
    ytempamp2 = ytempamp;
    ytempamp2(ytempamp<0) = 0;
    ytemptotamp = imgaussfilt( amplitudespectrums(:,:,[2:51,53:200])- mean(amplitudespectrums(:,:,[2:51,53:200]),3),2);
    yampthreshold = prctile(ytemptotamp(:),99.99);
  
    cmapbwr = buildcmap('rwb',64);
    figure(figurenum+6 ) 
    imagesc((xtempamp2-ytempamp2).*intensitymask) ; axis tight ;axis equal ;axis off 
    %colorbar ; set(gca,'FontSize',22)
    colormap(cmapbwr) ;caxis([-100 100]) ; %title([' x-y'])
    rectangle('Position',[1 1 255 255])


    
end


%% supfigure 2
% figure; 
% hold on
% boxplot(rawmatrix_off{1},'positions',1,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_off{3},'positions',2,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_off{5},'positions',3,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_off{7},'positions',4,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_off{9},'positions',5,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_on{1},'positions',1.4,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_on{3},'positions',2.4,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_on{5},'positions',3.4,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_on{7},'positions',4.4,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% boxplot(rawmatrix_on{9},'positions',5.4,'Whisker',0,'OutlierSize',3,'Jitter',0.2)
% label = {'2-off','2-on','4-off','4-on','6-off','6-on','8-off','8-on','10-off','10-on'};
% xlim([0.5 5.9])
% set(gca,'XTick',[1 1.4 2 2.4 3 3.4 4 4.4 5 5.4],'XTickLabel',label)

%% Violin Plot
% ON = Parallel; OFF = Perpendicular 
% Initialize the data values cell matrix 
data_vals = cell(10,1); 
% Initialize color 
clrs = cell(10,1); 
% Initialize labels (labels copied from above)
label = {'2-off','2-on','4-off','4-on','6-off','6-on','8-off','8-on','10-off','10-on'};
label_new = cell(size(label)); 
% Initialize counter 
cnt = 1;
for k = 1:10 
    % Even: Parallel i.e. on
    if mod(k,2) == 0 
        clrs{k,1} = 'c'; 
        data_vals{k,1} = rawmatrix_on{cnt};
        % Update the counter 
        cnt = cnt + 2; 
    % Odd: Perpendicular i.e. off
    else
        clrs{k,1} = 'm'; 
        data_vals{k,1} = rawmatrix_off{cnt};

    end 
    % Change the name of the label 
    newname = strrep(label{1,k},'off','\perp');
    newname = strrep(newname,'on','||'); 
    label_new{1,k} = newname; 
end 

% >>>> PLOT SETTINGS 
plot_settings = struct(); 
% Set the x-axis labels
plot_settings.xticklabel = label_new; 
% Change the rotation of the x-tick 
plot_settings.xtickrotation = -20; 
% Set the violin colors and transparency 
plot_settings.colorfill = clrs;
plot_settings.bordercolor = clrs; 
plot_settings.filltransparency = 1; 
% Line type for the median line 
plot_settings.linetype = {'-'};
% Change the y axis label 
plot_settings.ylabel = 'Magnitude (AU)'; 
% Change the x axis label 
plot_settings.xlabel = {'Distance away from bead center (\mu m)'...
    %,'and parallel or perpendicular to axis of oscillation'
    }; 
% Set the font sizes to be 12
plot_settings.font_size = 12; 
plot_settings.xlabelsize = 12; 
plot_settings.ylabelsize = 12; 
% Change width of the lines
plot_settings.borderwidth = 0.5; 
plot_settings.linewidth = 0.5;

% >>>> Generate the plots 
ax = figure; 
[cond_des, main_output, secondary_output, type] = ...
    plotViolinCell(data_vals, plot_settings); 
%% figure 8
% gradient
clear
close all 
%load('bead2-noosc.tif.mat')
%meanintimage = squeeze(mean(flashbackstack,3)) ;
%[x_gradient ,y_gradient ] =  gradient(meanintimage);

load('bead2-xosc.tif.mat')
meanintimage = squeeze(mean(flashbackstack,3)) ;
[x_gradient ,~ ] =  gradient(meanintimage);
load('bead2-yosc.tif.mat')
meanintimage = squeeze(mean(flashbackstack,3)) ;
[~,y_gradient ] =  gradient(meanintimage);



[cmap]=buildcmap('bwr',64);
figure(1)
imagesc(x_gradient); axis equal; axis tight; axis off; caxis([-1000 1000])
colormap(cmap) ; %colorbar; set(gca, 'FontSize',22);
hold on; rectangle('Position',[110,130,40,50]) ; hold off

figure(2)
imagesc(y_gradient); axis equal; axis tight; axis off; caxis([-1000 1000])
colormap(cmap) ; colorbar; set(gca, 'FontSize',22);
hold on; rectangle('Position',[110,130,40,50]) ; hold off

figure(3)
imagesc(x_gradient( 130:180,110:150)); axis equal; axis tight; axis off; caxis([-300 300])
colormap(cmap) ; colorbar; set(gca, 'FontSize',22);

figure(4)
imagesc(y_gradient( 130:180,110:150)); axis equal; axis tight; axis off; caxis([-300 300])
colormap(cmap) ; colorbar; set(gca, 'FontSize',22);

% load('bead2-xosc.tif.mat')
% tempamp= amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3);
% tempphase= angle(myffts(:,:,52));
%     
% figure(5)
% imagesc(tempamp( 130:180,110:150)); axis equal; axis tight; axis off; caxis([0 150])
% colormap jet ; %colorbar; set(gca, 'FontSize',22);
% 
% figure(6)
% imagesc(tempphase( 130:180,110:150)); axis equal; axis tight; axis off; %caxis([-300 300])
% colormap hsv ; %colorbar; set(gca, 'FontSize',22);

% load('bead2-yosc.tif.mat')
% tempamp= amplitudespectrums(:,:,52) - mean(amplitudespectrums(:,:,[2:51,53:200]),3);
% tempphase= angle(myffts(:,:,52));
%     
% figure(7)
% imagesc(tempamp( 130:180,110:150)); axis equal; axis tight; axis off; caxis([0 150])
% colormap jet ; colorbar; set(gca, 'FontSize',22);
% 
% figure(8)
% imagesc(tempphase( 130:180,110:150)); axis equal; axis tight; axis off; %caxis([-300 300])
% colormap hsv ; colorbar; set(gca, 'FontSize',22);



%% Supplemental videos 3,4,5,6
% manually cycle through the code depending on the comments
clear 
close all

for loop = 1:4

    if (loop==2||loop==4)
        load('bead2-xosc.tif.mat') % sup 4,6
    end
    if (loop==1||loop==3)
        load('bead2-noosc.tif.mat') % sup 3,5
    end
    if loop==2; filename = 'xosc.gif'; end %sup 4
    if loop==1; filename = 'noosc.gif';end %sup 3
    if loop==3; filename = 'noosc-phase.gif';end %sup 5
    if loop==4; filename = 'xosc-phase.gif'; end %sup 6

    h = figure(1);
    axis tight manual % this ensures that getframe() returns a consistent size
    for n = 1:200
        if (loop==3||loop==4)
            imagesc(angle(myffts(:,:,n))) %sup 5,6
        end
        if (loop==1||loop==2)
            imagesc(amplitudespectrums(:,:,n)) %sup 3,4
        end

        axis off equal
        cbar = colorbar;
             
        if (loop==3||loop==4)
            caxis([-pi pi]) %sup 5,6
            ylabel(cbar,'phase') %sup 5,6
            colormap hsv %sup 5,6
        end
        if (loop==1||loop==2)
            caxis([0 150]) %sup 3,4
            ylabel(cbar,'magnitude (AU)') %sup 3,4
        end
        
        title(['Frequency ' num2str(fq(n)) 'Hz'] )
        set(gca,'fontsize',12)
        truesize
        drawnow 
          % Capture the plot as an image 
          frame = getframe(h); 
          im = frame2im(frame); 
          [imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
          if n == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
    end

end


%% function for costom cmap
function [cmap]=buildcmap(colors,bitdepth)
% [cmap]=buildcmap(colors)
%
% This function can be used to build your own custom colormaps. Imagine if
% you want to display rainfall distribution map. You want a colormap which
% ideally brings rainfall in mind, which is not achiveved by colormaps such
% as winter, cool or jet and such. A gradient of white to blue will do the
% task, but you might also use a more complex gradient (such as
% white+blue+red or colors='wbr'). This function can be use to build any
% colormap using main colors rgbcmyk. In image processing, w (white) can be
% used as the first color so that in the output, the background (usually
% with 0 values) appears white. In the example of rainfall map, 'wb' will
% produce a rainfall density map where the background (if its DN values are
% 0) will appear as white.
%
% Inputs:
%  colors: string (char) of color codes, any sequence of rgbcmywk
%  representing different colors (such as 'b' for blue) is acceptable. If a
%  gradient of white to blue is needed, colors would be 'wb'; a rainbow of
%  white+blue+red+green would be 'wbrg'.
%
% Example:
%  [cmap]=buildcmap('wygbr');
% %try the output cmap:
% im=imread('cameraman.tif');
% imshow(im), colorbar
% colormap(cmap) %will use the output colormap
%
% First version: 14 Feb. 2013
% sohrabinia.m@gmail.com
%--------------------------------------------------------------------------
if nargin<1
    colors='wrgbcmyk';
end
if ~ischar(colors)
    error(['Error! colors must be a variable of type char with '...
        'color-names, such as ''r'', ''g'', etc., '...
        'type ''help buildcmap'' for more info']);
end
ncolors=length(colors)-1;
bins=round(bitdepth/ncolors);
% diff1=255-bins*ncolors;
vec=zeros(300,3);
switch colors(1)
    case 'w'
        vec(1,:)=1;
    case 'r'
        vec(1,:)=[1 0 0];
    case 'g'
        vec(1,:)=[0 1 0];
    case 'b'
        vec(1,:)=[0 0 1];
    case 'c'
        vec(1,:)=[0 1 1];
    case 'm'
        vec(1,:)=[1 0 1];
    case 'y'
        vec(1,:)=[1 1 0];
    case 'k'
        vec(1,:)=[0 0 0];
end
for i=1:ncolors
 beG=(i-1)*bins+1;
 enD=i*bins+1; %beG,enD
 switch colors(i+1)
     case 'w'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD,
     case 'r'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
     case 'g'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
     case 'b'         
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
     case 'c'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
     case 'm'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
     case 'y'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
     case 'k'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
 end
end
cmap=vec(1:bins*ncolors,:);
end %end of buildcmap

%% Function for calculating power spectrum. 
function [fq,amplitudespectrum,myfft] = myamplitudespectrum(data,frequency)
    newlength=length(data);
    myfft = fft(data);
    amplitudespectrum= abs(myfft./newlength);
    freqHz = (0:1:newlength-1).*frequency./newlength;
    fq = freqHz(1:floor(newlength/2));
    amplitudespectrum= amplitudespectrum(1:floor(newlength/2));
    amplitudespectrum(2:end)= amplitudespectrum(2:end)*2; %single sided
%figure
%loglog(freqHz(1:length(data)/2), amplitudespectrum(1:length(data)/2))
end
