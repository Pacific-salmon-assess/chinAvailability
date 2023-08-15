===================================================================
SRKW FORAGING DATA PACKAGE FOR CAM FRESHWATER (17-Nov-2022)
===================================================================

Data included in this package:
==============================
(1) Rasters providing the probability of CFA (common foraging area)
	> Common foraging area = where the probability of foraging exceeds 25%
	> One raster provided for each region (Haro & Swiftsure)
	> Raster datum/projection: NAD83 BC Albers
	> Filenames: Swiftsure_post_forage_exc_0.25_NAD83_BCAlbers.tif, Haro_post_forage_exc_0.25_NAD83_BCAlbers.tif
(2) Rasters providing the probability of FFA (frequent foraging area)
	> Frequent foraging area = where the probability of foraging exceeds 50%
	> One raster provided for each region (Haro & Swiftsure)
	> Raster datum/projection: NAD83 BC Albers
	> Filenames: Swiftsure_post_forage_exc_0.5_NAD83_BCAlbers.tif, Haro_post_forage_exc_0.5_NAD83_BCAlbers.tif
(3) Shapefiles delineating CFAs of high confidence
	> Where probability of foraging exceeds 25% in more than 70/80/90% of model iterations
	> 3 shapefiles provided for each region (Haro & Swiftsure), for 70, 80, and 90% confidence bounds, named "0.7prop", 0.8prop", and "0.9prop" respectively.
	> Filenames:	haro.forage.0.25exc.0.7prop.poly_NAD83_BCAlbers.shp,
			haro.forage.0.25exc.0.8prop.poly_NAD83_BCAlbers.shp,
			haro.forage.0.25exc.0.9prop.poly_NAD83_BCAlbers.shp,
			swiftsure.forage.0.25exc.0.7prop.poly_NAD83_BCAlbers.shp,
			swiftsure.forage.0.25exc.0.8prop.poly_NAD83_BCAlbers.shp,
			swiftsure.forage.0.25exc.0.9prop.poly_NAD83_BCAlbers.shp
(4) Shapefiles delineating FFAs of high confidence
	> Where probability of foraging exceeds 50% in more than 70/80/90% of model iterations
	> 3 shapefiles provided for each region (Haro & Swiftsure), for 70, 80, and 90% confidence bounds, named "0.7prop", 0.8prop", and "0.9prop", respectively.
	> Filenames:	haro.forage.0.5exc.0.7prop.poly_NAD83_BCAlbers.shp,
			haro.forage.0.5exc.0.8prop.poly_NAD83_BCAlbers.shp,
			haro.forage.0.5exc.0.9prop.poly_NAD83_BCAlbers.shp,
			swiftsure.forage.0.5exc.0.7prop.poly_NAD83_BCAlbers.shp,
			swiftsure.forage.0.5exc.0.8prop.poly_NAD83_BCAlbers.shp,
			swiftsure.forage.0.5exc.0.9prop.poly_NAD83_BCAlbers.shp

SEE FIGURES INCLUDED PACKAGE visualizing the CFA and FFA rasters with polygon shapefiles overlaid:
> SRKWBehav_Fig5_common foraging areas_labs.png
> SRKWBehav_Fig6_frequent foraging areas_labs.png
