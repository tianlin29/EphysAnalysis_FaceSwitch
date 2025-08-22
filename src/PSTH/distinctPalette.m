
function palette = distinctPalette(index)
% 
% 
% 
%     color_lib = distinctPalette;
%     figure;
%     hold on;
%     for i = 1 : 16
%         plot(1:100,(1:100)+i*5, 'Color', color_lib{i}, 'LineWidth', 1);
%     end
% 

% 
% RK, 1/4/2017
% 


color_lib = {
    [240,163,255]/255
    [0,117,220]/255
    [153,63,0]/255
    [0,0,0]/255
    [0,92,49]/255
    [43,206,72]/255
    [128,128,128]/255
    [143,124,0]/255
    [157,204,0]/255
    [194,0,136]/255
    [0,51,128]/255
    [255,164,5]/255
    [66,102,0]/255
    [255,0,16]/255
    [0,153,143]/255
    [116,10,255]/255
};

if nargin<1 || isempty(index)
    palette = color_lib;
else
    palette = color_lib(index);
end


