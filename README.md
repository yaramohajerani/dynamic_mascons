Dynamic Mascons
====================

[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/yaramohajerani/FrontLearning/blob/master/LICENSE)
[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/yaramohajerani/dynamic_mascons/blob/main/spherical_voronoi.ipynb)
[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/yaramohajerani/dynamic_mascons/blob/main/voronoi_to_mascon.ipynb)

### By Yara Mohajerani

## Summary 

This respository showcases global dynamic mascon configurations that are regionally optimizied through a series of designed fixed-points and an iterative re-tessellation scheme based on polygon centroids. 

In the first notebook, `spherical_voronoi.ipynb`, I go through the steps of constructing a spherical voronoi tessellation with a set of fixed points. I choose two points on the Karakoram region in High Mountain Asia as our fixed points. 
[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/yaramohajerani/dynamic_mascons/blob/main/spherical_voronoi.ipynb)

In the second notebook, the resulting final configuration is used to construct a harmonic mascon representation, which is then used to calculate the sensitivity kernel of the mascons resulting from the inversion (refer to [Jacob et al. 2011](http://doi.org/10.1007/s00190-011-0522-7)). The resulting kernel for the fixed points is localized and minimal leakage. This is in contrast to many configurations with regional coverage that are prone to divergence of the kernel at the boundaries of the domain, which illustrates the effectiveness of this method.
[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/yaramohajerani/dynamic_mascons/blob/main/voronoi_to_mascon.ipynb)

## More Background

Some of my past work has focused on the development of regionally-optimizied GRACE mascons. E.g.
* [Optimizied Mascons for Totten glacier in Antarctica](https://doi.org/10.1029/2018GL078173)
* [Optimzied Mascons for Getz and Amery basins in Antarctica](https://doi.org/10.1029/2019GL084665)

However, there is a need for a more robust approach to create *global* mascon configurations that are dynamic and depend on the geophysical feature of the region of interest to focus on. In this repo we explore and develop the framework to create dynamic mascon configurations based an underlying desired density.

This draws on some previous work based on an iterative Voronoi Tessallation approach that I developed for Antarctica:

![](./imgs/sample_antarctica.gif)

However, the extension this to a global definition requires a 3D geodetic grid with no distortion.

Some useful resources:
* [Hexagonal CNN](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8853238)
* [HexASCII](https://onlinelibrary.wiley.com/doi/epdf/10.1111/tgis.12304)
* [SciPy spherical voronoi module](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.SphericalVoronoi.html)