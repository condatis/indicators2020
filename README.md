# Testing connectivity indicators in R

[www.github.com/condatis/indicators2020](www.github.com/condatis/indicators2020), Licensed with Open Government License v.3 (OGL), see License.txt

Jenny A. Hodgson, 2020

Contact: jenny.hodgson\@liverpool.ac.uk; github\@condatis.org.uk

This repository contains a number of R functions which facilitate the calculation of connectivity indicators for a region and habitat type of interest. Condatis 'speed' (see [www.condatis.org.uk](www.condatis.org.uk)) is one of the indicators; two moving-window-based indicators are also included. In order to test the behaviour of the indicators, functions are included that increase the proportion of habitat in a raster landscape, using a number of different weighting schemes to achieve different spatial patterns.

The code was written with the intention of contributing to a collaborative project with the UK Centre for Ecology and Hydrology, funded by Defra (the UK department of environment, food and rural affairs). The collaboration led to a paper [^1] and two reports to Defra. Please note that, although much of this code has potential to be re-used in projects other than the one for which the scripts were intended, it was not written with broad re-usability in mind. Only minimal comments are included for guidance. This repository serves to differentiate functions that were authored by me, from others authored collaboratively, and edited after 2020, and therefore it does not contain everything necessary to reproduce the products that UKCEH delivered to Defra. Anyone wishing to use the code is encouraged to contact me to understand it better.

Preferred attribution if you re-use this work under the OGL: "derived from" Hodgson (2020) 'Testing connectivity indicators in R' www.github.com/condatis/indicators2020, DOI: 10.5281/zenodo.7950780

[^1]:Mancini, Francesca, Jenny A. Hodgson, and Nick J. B. Isaac. 2022.“Co-Designing an Indicator of Habitat Connectivity forEngland.” Frontiers in Ecology and Evolution 10 (July).https://doi.org/10.3389/fevo.2022.892987.
## Rapid reference - list of functions

In GIS_functions.R

`prepregion()`

`makehabitat()`

`makehabitatpoints()`

In neighbwindowfuns2.R

`hanski()`

`expweightmat()`

`ccvc()`

In condatisfunctions.R

`st4directions()`

`speednflow()`

In patchiness.R

`patchisolation()`

`thumb1()`

In addinghabitat.R

`addhabitat3ways()`

In addinghabitat_lowish.R

`addhabitat.lowish()`
