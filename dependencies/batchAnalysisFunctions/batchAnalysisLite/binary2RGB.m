function imageMat = binary2RGB(logicMat, rgbVal);

% if max(rgbVal) <= 1
%   rgbVal = rgbVal * 255;
% end

% Cycle through the layers, adding the rgb value where 1s are.
imageMat = ones([size(logicMat),3]);
for layer_i = 1:3
  imageMat(:,:,layer_i) = logicMat * rgbVal(layer_i);
end

end