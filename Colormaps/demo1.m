%% Colormaps that are pretty awesome: DEMO1
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
close all;

screen=get(0,'ScreenSize') ;
w0=screen(1);
h0=screen(2);
w =screen(3);
h =screen(4);
%% Generate sample data
load flujet
X=X.';

%% Plot original with jet and parula
% Parula
h2=figure('name','Sample data with "parula" colormap');
set(h2,'position',[w0,h0,w/3,h/2])
imagesc(X)

if verLessThan('matlab', '8.4')
   % if parula is the "future"
   colormap(parula());
else
   % if parula is already in Matlab
   colormap('parula');
end
axis image
xlabel('Parula Colormap')
set(gca,'xtick',[],'ytick',[]) % This is axis off without offing the labels

% Jet
h1=figure('name','Sample data with "jet" colormap');
set(h1,'position',[w0,h/2,w/3,h/2])
imagesc(X)
colormap('jet')
xlabel('jet Colormap')
set(gca,'xtick',[],'ytick',[]) % This is axis off without offing the labels
axis image


%% Load new colormaps

m=40
cmA=py_A_cmap(m);
cmB=py_B_cmap(m);
cmC=py_C_cmap(m);
cmD=viridis(m);


%% Plot new colormaps
h3=figure('name','Super-cool new colormaps that you can easily use');
set(h3,'position',[w/3,h0,2*w/3,h])
if verLessThan('matlab', '8.4')
% If you are using old Matlab figure engine do  it this way
% (some very lousy colormap problems before)
    subplot(2,2,1,'Position',[0.05 0.55 0.4 0.4])
    subimage(uint8(X/max(X(:))*255),cmA)
    xlabel('python A Colormap')
    set(gca,'xtick',[],'ytick',[]) 


    subplot(2,2,2,'Position',[0.55 0.55 0.4 0.4])
    subimage(uint8(X/max(X(:))*255),cmB)
    xlabel('python B Colormap')
    set(gca,'xtick',[],'ytick',[]) 


    subplot(2,2,3,'Position',[0.05 0.05 0.4 0.4])
    subimage(uint8(X/max(X(:))*255),cmC)
    xlabel('python C Colormap')
    set(gca,'xtick',[],'ytick',[])
    
    subplot(2,2,4,'Position',[0.55 0.05 0.4 0.4])
    subimage(uint8(X/max(X(:))*255),cmD)
    xlabel('python D Colormap')
    set(gca,'xtick',[],'ytick',[])
else
    sp1=subplot(2,2,1,'Position',[0.05 0.55 0.4 0.4]);
    imagesc(X)
    colormap(sp1,cmA)
    xlabel('python A Colormap')
    set(gca,'xtick',[],'ytick',[])
    
    sp2=subplot(2,2,2,'Position',[0.55 0.55 0.4 0.4]);
    imagesc(X)
    colormap(sp2,cmB)
    xlabel('python B Colormap')
    set(gca,'xtick',[],'ytick',[])

    sp3=subplot(2,2,3,'Position',[0.05 0.05 0.4 0.4]);
    imagesc(X)
    colormap(sp3,cmC)
    xlabel('python C Colormap')
    set(gca,'xtick',[],'ytick',[])
    
    sp4=subplot(2,2,4,'Position',[0.55 0.05 0.4 0.4]);
    imagesc(X)
    colormap(sp4,cmD)
    xlabel('python D Colormap')
    set(gca,'xtick',[],'ytick',[])
end
