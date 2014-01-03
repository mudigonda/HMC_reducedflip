function [ dEdX ] = dEdX_rough( X, scale1, scale2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sinX = sin(X*2*pi/scale2);
dEdX = X./scale1^2 + -sinX*2*pi./scale2;

end

