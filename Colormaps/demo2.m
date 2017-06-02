%% Colormaps that are pretty awesome: DEMO2
%  In this demo we will show how to supereasily use the new colormaps
%
%
% this is an exampe of 4 colormaps that are being considered as the default
% colormap in python's matplotlib lybrary.
%
% All of them look quite good and they dont have any official name, so at
% the moment they are A,B,C,D.
%
% colormaps from https://github.com/bids/colormap
%
% Ander Biguri
%% Clear workspace and get screen data
clear;
clc
% close all;


%% Generate sample data
X=peaks(200);

%% Load Colomaps

jet=colormap('jet');
parula=parula();
pyA=py_A_cmap();
pyB=py_B_cmap();
pyC=py_C_cmap();
pyD=viridis();


%% Plot

surf(X,'linestyle','none');
%% Chose colormap
% Use 1 only, else it will  just use the last
% CTRL+R -> comment line
% CTRL+T -> Uncomment line

% colormap(jet);
% colormap(parula);
% colormap(pyA);
%  colormap(pyB);
%  colormap(pyC);
 colormap(pyD);

