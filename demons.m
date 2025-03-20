% demons.m
%
% Purpose: Create transformed version of input image A to closest match
% with input image B using image registration (phase correlation, demons)
% 
% Author(s): Jon Seibert
% Last updated: 13 Nov 2023
% 
% Inputs: MxN Binary array,  MxN lon array, MxN lat array
% Outputs: MxN NP array
% 
% Dependencies: Matlab image processing toolbox
%%

function a_trans = demons(data_a,data_b)

tformEstimate = imregcorr(data_a,data_b,"translation");
Rfixed = imref2d(size(data_b));
a_corr = imwarp(data_a,tformEstimate,"OutputView",Rfixed);

[D,a_trans] = imregdemons(a_corr,data_b);

%a_trans(a_trans == 0) = -30;

%%

end
