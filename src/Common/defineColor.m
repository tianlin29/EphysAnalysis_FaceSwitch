function [redblue, clims] = defineColor()
m = 256;
redblue = [linspace(0,1,m/2)', linspace(0,1,m/2)', ones(m/2,1); ...
    ones(m/2,1), linspace(1,0,m/2)', linspace(1,0,m/2)'];
clims = [-1 1];
end
