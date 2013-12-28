function [ E ] = E_rough( X, scale1, scale2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cosX = cos(X*2*pi/scale2);
E = sum( (X.^2) / (2*scale1^2) + cosX );

end

