function [ score ] = fScore( sample, template )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sample=mat2gray(sample);
template=mat2gray(template);

commonEl=sample.*template;
precision=sum(commonEl(:))/sum(sample(:));
recall=sum(commonEl(:))/sum(template(:));
score=2*precision*recall/(precision+recall);
end

