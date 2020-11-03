Dynamic Mascons
====================

[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/yaramohajerani/FrontLearning/blob/master/LICENSE)
[![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/yaramohajerani/dynamic_mascons/blob/main/spherical_voronoi.ipynb)

### By Yara Mohajerani

Some of my past work has focused on the development of regionally-optimizied GRACE mascons. E.g.
* [Optimizied Mascons for Totten glacier in Antarctica](https://doi.org/10.1029/2018GL078173)
* [Optimzied Mascons for Getz and Amery basins in Antarctica](https://doi.org/10.1029/2019GL084665)

However, there is a need for a more robust approach to create *global* mascon configurations that are dynamic and depend on the geophysical feature of the region of interest to focus on. In this repo we explore and develop the framework to create dynamic mascon configurations based an underlying desired density.

This draws on some previous work based on an iterative Voronoi Tessallation approach that I developed for Antarctica:

![](./imgs/sample_antarctica.gif)

However, the extension this to a global definition requires a 3D geodetic grid with no distortion, e.g. ![](https://i.stack.imgur.com/L6kmP.jpg)

Some useful resources:
* [Hexagonal CNN](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8853238)
* [HexASCII](https://onlinelibrary.wiley.com/doi/epdf/10.1111/tgis.12304)
* [SciPy spherical voronoi module](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.SphericalVoronoi.html)
