function [ D ] = imreconstructclosing( m,g )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

m = imcomplement(m);
g = imcomplement(g);

D = imreconstruct(m,g);

D = imcomplement(D);
end

