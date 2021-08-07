# Automated Dynamic Mascon Generation for GRACE and GRACE-FO Harmonic Processing

This repository presents a novel automated methodology for creating dynamic global mascons using iterative spherical Voronoi tessellations for the processing of GRACE (Gravity Recovery and Climate Experiment) and GRACE-FO (Follow-On) harmonics.

Click the image below for an interactive dynamic global configurations for the Karakoram region in High Mountain Asia:

<p><a href="karakoram_config.html">
<img src="tessellations_karakoram.png" alt="Karakoram Configuration">
</a></p>

The associated publication can be found at

> Mohajerani, Y.; Shean, D.; Arendt, A.; Sutterley, T.C. Automated Dynamic Mascon Generation for GRACE and GRACE-FO Harmonic Processing. Remote Sens. 2021, 13, 3134. [https://doi.org/10.3390/rs13163134](https://doi.org/10.3390/rs13163134)

### Abstract

> Commonly used mass-concentration (mascon) solutions estimated from Level-1B Gravity Recovery and Climate Experiment (GRACE) and GRACE Follow-On data, provided by processing centers such as the Jet Propulsion Laboratory (JPL) or the Goddard Space Flight Center (GSFC), do not give users control over the placement of mascons or inversion assumptions, such as regularization. While a few studies have focused on regional or global mascon optimization from spherical harmonics data, a global optimization based on the geometry of geophysical signal as a standardized product with user-defined points has not been addressed. Finding the optimal configuration with enough coverage to account for far-field leakage is not a trivial task and is often approached in an ad-hoc manner, if at all. Here, we present an automated approach to defining non-uniform, global mascon solutions that focus on a region of interest specified by the user, while maintaining few global degrees of freedom to minimize noise and leakage. We showcase our approach in High Mountain Asia (HMA) and Alaska, and compare the results with global uniform mascon solutions from range-rate data. We show that the custom mascon solutions can lead to improved regional trends due to a more careful sampling of geophysically distinct regions. In addition, the custom mascon solutions exhibit different seasonal variation compared to the regularized solutions. Our open-source pipeline will allow the community to quickly and efficiently develop optimized global mascon solutions for an arbitrary point or polygon anywhere on the surface of the Earth.


The Github repo can be cited as
> Mohajerani, Y. Yaramohajerani/Dynamic_mascons: Release for Accepted Manuscript (v2.0). 2021. Available online: 
[https://doi.org/10.5281/zenodo.5167967](https://doi.org/10.5281/zenodo.5167967)



#### License
The code presented in this repository is under the [MIT license](./../LICENSE).
Copyright &copy; 2021 Yara Mohajerani
