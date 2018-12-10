# EMODnet-Biology-Zooplankton-Baltic

#  Input data

* Zookplankton observations for 40 different species from the Baltic from the Swedish SHARK database from EMODnet Biology, the Finnish data from the [NOAA Copepod database](https://www.st.nmfs.noaa.gov/copepod/), the German and Polish from the HELCOM DOME database. By the finishing date of the project, all these dataset will be integrated in the EMODnet Biology database.
* [Dissolved oxygen](http://www.emodnet-chemistry.eu/products/catalogue#/metadata/087a72c0-c243-11e8-bac2-5ce0c5469bc7) (from EMODnet Chemistry)
* [Salinity](https://www.seadatanet.org/Products#/metadata/bf35a7c5-c843-4a23-8040-07ddcf3d8e71) (from SeaDataNet)
* [Temperature](https://www.seadatanet.org/Products#/metadata/bf35a7c5-c843-4a23-8040-07ddcf3d8e71) (from SeaDataNet)
* [https://oceancolor.gsfc.nasa.gov/l3/](Chlorophyll concentration) (MODIS-Aqua from NASA)
* [Bathymetry](https://www.gebco.net/) (from GEBCO)
* [Distance from coast](https://gcmd.nasa.gov/KeywordSearch/Metadata.do?Portal=idn_ceos&KeywordPath=%5BData_Center%3A+Short_Name%3D%27PacIOOS%27%5D&OrigMetadataNode=GCMD&EntryId=dist2coast_1deg&MetadataView=Full&MetadataType=0&lbnode=mdlb1) (from GSFC, NASA)


# Data product description

A gridded data product for 40 zooplankton species using DIVAnd[1,2] and the neural network library Knet[3]. The neural network uses the variables dissolved oxygen, salinity, temperature, chlorophyll concentration bathymetry and the distance from coast as input. Additionally the position (latitude and longitude) and the year are provided to the neural network.

Abundance values in the NetCDF files are expressed in number per m² and transformed by the function log(x/a + 1) where a is 1 m⁻².

The covers the area from  9°E to 30.8°E and 53°N to 66.1°N at a resolution of a tenth of a degree. Gridded data product for the years 2007, 2008, 2010, 2011, 2012 and 2013 have been made. No observations were available for the year 2009. The fields represent the yearly average abundance.

For every specie the correlation length and signal to noise ratio is estimated using the spatial variability of the observations.

The interpolation error is computed using the so-called "clever poor man's error" [4].


The full list of the species is:
* Acartia (Acanthacartia) bifilosa
* Acartia (Acanthacartia) tonsa
* Acartia (Acartiura) clausi
* Acartia (Acartiura) longiremis
* Amphibalanus improvisus
* Appendicularia
* Bivalvia
* Bosmina (Eubosmina) coregoni
* Bryozoa
* Calanus finmarchicus
* Centropages
* Centropages hamatus
* Cercopagis (Cercopagis) pengoi
* Cnidaria
* Cyclopoida
* Daphnia
* Daphnia cristata
* Echinodermata
* Eurytemora
* Evadne nordmanni
* Fritillaria borealis
* Gastropoda
* Harpacticoida
* Keratella cochlearis
* Keratella cruciformis
* Keratella eichwaldi
* Keratella quadrata
* Limnocalanus macrurus macrurus
* Mysidae
* Oikopleura
* Oithona
* Paracalanus
* Pleopis polyphemoides
* Podon intermedius
* Podon leuckartii
* Polychaeta
* Pseudocalanus
* Rotifera
* Synchaeta
* Temora longicornis


## References

[1] A. Barth, J.-M. Beckers, C. Troupin, A. Alvera-Azcárate, and L. Vandenbulcke. divand-1.0: n-dimensional variational data analysis for ocean observations. Geoscientific Model Development, 7(1):225–241, 2014. doi: [10.5194/gmd-7-225-2014](https://doi.org/10.5194/gmd-7-225-2014).

[2] C. Troupin, A. Barth, D. Sirjacobs, M. Ouberdous, J.-M. Brankart, P. Brasseur, M. Rixen, A. Alvera-Azcárate, M. Belounis, A. Capet, F. Lenartz, M.-E. Toussaint, and J.-M. Beckers. Generation of analysis and consistent error fields using the Data Interpolating Variational Analysis (DIVA). Ocean Modelling, 52–53:90–101, 2012. doi: [10.1016/j.ocemod.2012.05.002](https://doi.org/10.1016/j.ocemod.2012.05.002).

[3] Yuret, D. Knet: beginning deep learning with 100 lines of julia. In Machine Learning Systems Workshop at NIPS 2016.

[4] J.-M. Beckers, A. Barth, C. Troupin, and A. Alvera-Azcárate. Approximate and efficient methods to assess error fields in spatial gridding with data interpolating variational analysis (DIVA). Journal of Atmospheric and Oceanic Technology, 31:515–530, 2014. doi: [10.1175/JTECH-D-13-00130.1](https://doi.org/10.1175/JTECH-D-13-00130.1).

# Code

The data products was made using DIVAnd in version 2.1.1 (https://github.com/gher-ulg/DIVAnd.jl) and Knet in version 0.9.2 (https://github.com/denizyuret/Knet.jl)
and Julia 0.6.4. See the REQUIRE file for a fill list of dependencies. The main entry point of this code is the file `emodnet_bio3.jl`. All input data is assumed to be in the directory `$HOME/tmp/Emodnet-Bio` (or the directory defined in the environement variable `DATADIR`)

<!--
Link to "Baltic Zooplankton Data Preparation" from Peter M.J. Herman and Lisa Sundqvist
# Link to data products

* within EMODnet geoviewer
* As netCDF: https://dox.ulg.ac.be/index.php/s/GTjQHky8I5zSHLF/download (to be transferred to VLIZ server)
-->
