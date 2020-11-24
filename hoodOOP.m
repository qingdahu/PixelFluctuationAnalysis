% hoodOOP - calculate the orientational order parameter in the neighborhood
% surrounding each pixel in an image 
%
%
% Arguments:
%   I       - orientation angles in radians with no NaN values
%               Class Support: numeric 
%   sze     - number of pixels in the horizontal and vertical directions to
%               to include. The area of the neighborhood is (2*sze+1)^2 
%     
%
% 
% Returns:
%   Ioop    - orientational order parameter for each pixel and its
%               neighbors
%               Class Support: double 
%   Ineigh  - number of pixels in the neighborhood of each pixel 
%               Class Support: double  
%
% Dependencies: 
%   calculate_OOP.m
%   MATLAB Version >= 9.5 
%
%
%   Written by: Tessa Morris
%   Advisor: Anna (Anya) Grosberg, Department of Biomedical Engineering 
%   Cardiovascular Modeling Laboratory 
%   University of California, Irvine 

function [Ioop, Ineigh] = hoodOOP(I,sze)
% Get the size of the image 
[max_d1,max_d2] = size(I); 

% Save the OOP image 
Ioop = zeros(size(I)); 
% Save the number of neighboring pixels accounted for 
Ineigh = zeros(size(I)); 

% Set boolean dontRun to false 
dontRun = false; 

% Make sure that size is a positive integer 
if mod(sze,1) ~= 0 
    sze = round(sze); 
    disp('The neighborhood size was rounded to an integer.'); 
end 
if sze < 1
    disp('The neighborhood size must be positive. No analysis possible'); 
    dontRun = true; 
end 

% Loop through all of the pixels in the image 
if ~dontRun 
    for poi_d1 = 1:max_d1
        for poi_d2 = 1:max_d2 
            % Get the starting position 
            start_d1 = poi_d1 - sze; 
            start_d2 = poi_d2 - sze;
            % Get the stopping positions 
            stop_d1 = poi_d1 + sze; 
            stop_d2 = poi_d2 + sze;

            % Correct the image section dimensions 
            if start_d1 < 1
                start_d1 = 1; 
            end 
            if start_d2 < 1
                start_d2 = 1; 
            end 
            if stop_d1 > max_d1
                stop_d1 = max_d1; 
            end 
            if stop_d2 > max_d2
                stop_d2 = max_d2; 
            end 

            % Get the temporary image 
            temp_im = I(start_d1:stop_d1, start_d2:stop_d2);
            disp(size(temp_im));
            % Remove any NaN
            temp_im = temp_im(:); 
            temp_im(isnan(temp_im)) = []; 
            % Calculate the OOP of the temporary image
            [ OOP, ~, ~, ~ ] = calculate_OOP( temp_im );
            Ioop(poi_d1,poi_d2) = OOP;    
            
            % Save the number of pixels that was used to calculate the oop 
            Ineigh(poi_d1,poi_d2) = length(temp_im) - 1;
        end 
    end
end

end

