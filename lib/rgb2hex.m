%======================================================================
%>@brief Converts RGB color triplets to hexadecimal format
%>
%>@param rgb (@b array) RGB color triplet
%>@retval hex (@b char) Hexadecimal color code
%>
%======================================================================
function [ hex ] = rgb2hex(rgb)
%% Check inputs: 

assert(nargin==1,'This function requires an RGB input.') 
assert(isnumeric(rgb)==1,'Function input must be numeric.') 

sizergb = size(rgb); 
assert(sizergb(2)==3,'rgb value must have three components in the form [r g b].')
assert(max(rgb(:))<=255& min(rgb(:))>=0,'rgb values must be on a scale of 0 to 1 or 0 to 255')

%% If no value in RGB exceeds unity, scale from 0 to 255: 
if max(rgb(:))<=1
    rgb = round(rgb*255); 
else
    rgb = round(rgb); 
end

%% Convert (Thanks to Stephen Cobeldick for this clever, efficient solution):

hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).'; 
hex(:,1) = '#';


end
%======================================================================
%>@file rgb2hex.m
%>@brief Conversion of color triplets to hexadecimal format
%>@note
%>External function, available at MATLAB Central.
%>@author Chad A. Greene
%>@date April of 2014
%>@see
%>MATLAB central<br>
%>[rgb2hex and hex2rgb](https://www.mathworks.com/matlabcentral/fileexchange/46289-rgb2hex-and-hex2rgb)
%>
%======================================================================