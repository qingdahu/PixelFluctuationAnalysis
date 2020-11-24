% May 15, 2020
% Plase make sure you select all the .TIF files you are going to need for
% the analysis.
clear all
close all
% select all .tif (do not select the linescans)
[filenames,location] = uigetfile('.tif','MultiSelect','on'); 
for filenum = 1:length(filenames)
    filename = char(filenames(filenum))  ;     
    flashbackstack = readinstack(filename,location); % Runs function below
    totalrows = size(flashbackstack,1); %should be 256
    totalcols = size(flashbackstack,2); %should be 256
    totalframes = size(flashbackstack,3); %should be 400
    frametime = 0.064904; % In seconds
    for row =1:totalrows
        for col =1:totalcols
            [fq,amplitudespectrum,myfft] = myamplitudespectrum(flashbackstack(row,col,:),1/frametime); %runs the function below
            amplitudespectrums(row,col,:) = amplitudespectrum;        
            myffts(row,col,:) = myfft;
        end
    end
    save([filename '.mat']) 
end

%% Function for calculating power spectrum. 
function [fq,amplitudespectrum,myfft] = myamplitudespectrum(data,frequency)
    newlength=length(data);
    myfft = fft(data);
    amplitudespectrum= abs(myfft./newlength);
    freqHz = (0:1:newlength-1).*frequency./newlength;
    fq = freqHz(1:floor(newlength/2));
    amplitudespectrum= amplitudespectrum(1:floor(newlength/2));
    amplitudespectrum(2:end)= amplitudespectrum(2:end)*2; %single sided
    %figure; loglog(freqHz(1:length(data)/2), amplitudespectrum(1:length(data)/2))
end

%% Function for reading in tiff stacks
function imagestack=readinstack(file,path)
    tiff_info = imfinfo([path,file]); % return tiff structure, one element per image
    imagestack = imread([path,file], 1) ; % read in first image
    %concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_info, 1)
        temp_tiff = imread([path,file], ii); %read in the next image
        imagestack = cat(3 , imagestack, temp_tiff); %save to the stack in matlab
    end
end