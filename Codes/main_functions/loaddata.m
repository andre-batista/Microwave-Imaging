function [model,es,ei,gs,gd] = loaddata(filepath,filename)
    load([filepath,filename,'_data.mat'])