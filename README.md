# xATLy: E*x*tending ICESat-2 **A**TL03 and ATL08 Products for **L**and and Vegetation Anal**y**ses
## Overview
The xATLy toolbox facilitates the reading, visualization, and processing of ICESat-2 ATL03 and ATL08 in HDF5 format with the end goal of generating customized along-track terrain and canopy height data sets.
<p align="center">
<img src = "https://github.com/Oht0nger/xATLy/blob/main/code/doc/xatly_ch.png" width = "900" />
</p>

## System Requirements and toolbox installation
* The toolbox requires MATLAB (Release 2021+), the Curve fitting toolbox, Mapping toolbox and the Statistics toolbox
* To install the toolbox, download the [package file](https://github.com/Oht0nger/xATLy/blob/main/toolbox/xATLy.mltbx) and follow [MATLAB'S documentation](https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html)
## Features
Overview of the toolbox features include:
* Read and write raw ATL03 point cloud
* Derive classified ATL03 point cloud 
* Derive aboveground level (normalized) point clouds
* Read ATL08 segment-level, canopy and terrain attributes
* Plot ATL03 point clouds
* Calculate custom canopy and terrain heights
* Generate and write auxillary geospatial layers as shapefile or KML
  * Reference ground tracks
  * 100-m ATL08 segment-level polygons
  * Custom segment-level polygons
## Examples
Explore the following examples on how to use the toolbox features.
* Working with ICE2Veg toolbox | [MATLAB Live Script](https://github.com/Oht0nger/xATLy/blob/main/code/doc/Working%20with%20xATLy%20Toolbox.mlx) | [pdf](https://github.com/Oht0nger/xATLy/blob/main/code/doc/Working%20with%20xATLy%20Toolbox.pdf)
* Custom Canopy and Terrain Height Calculations | [MATLAB Live Script](https://github.com/Oht0nger/xATLy/blob/main/code/doc/Custom%20Canopy%20and%20Terrain%20Height%20Calculation.mlx) | [pdf](https://github.com/Oht0nger/xATLy/blob/main/code/doc/Custom%20Canopy%20and%20Terrain%20Height%20Calculation.pdf)
