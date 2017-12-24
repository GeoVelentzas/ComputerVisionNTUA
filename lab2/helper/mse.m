function [ error ] = mse( sample,template )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
sample=mat2gray(sample);
template=mat2gray(template);
A=~sample.*template;
B=sample.*~template;
error=0.5*(sum(A(:))/sum(template(:))+sum(B(:))/sum(template(:)));

end

