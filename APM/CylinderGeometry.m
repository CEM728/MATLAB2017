classdef CylinderGeometry
  properties
    % simulation properties
    orders@double;

    % background properties
    ep@double scalar complex;
    mu@double scalar complex;

    % field properties
    k@double scalar complex;
    beta@double scalar complex;

    % cylinder properties
    a@double scalar;
    epi@double scalar complex;
    mui@double scalar complex;
  end
end
