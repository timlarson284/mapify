# Data Format Control dictionary
# this is mainly used for holding the confidence values and such used in the
# land cover layers.
# lc -> Land Cover layers
# lcc -> Land Cover Confidence layers
# lccf -> Land Cover Confidence layers (fill techniques)
dfc = {'lc_insuff': 10,  # insufficient data, such as at the end of a time series
       'lcc_growth': 151,
       'lcc_decline': 152,
       'lccf_nomodel': 201,
       'lccf_forwards': 202,
       'lccf_samelc': 211,
       'lccf_difflc': 212,
       'lccf_back': 213,
       'lccf_afterbr': 214,
       }

chg_magbands = ('green', 'red', 'nir', 'swir1', 'swir2')
chg_begining = '1982-01-01'

# ts_start: '1982-01-01'
# ts_stop: '2018-01-01'
#
# nodata: 0
# disturbed: 9
# series-end: 10
#
# conus-tileaff: [-2565585, 150000, 0 , 3314805, 0, -150000]
# conus-chipaff: [-2565585, 3000, 0 , 3314805, 0, -3000]
#
# conus-extent:
#   xmin: -2565585
#   ymax: 3314805
#   xmax: 2384415
#   ymin: 14805
# alaska-extent:
#   xmin: -851715
#   ymax: 2474325
#   xmax: 1698285
#   ymin: 374325
# hawaii-extent:
#   xmin: -444345
#   ymax: 2168895
#   xmax: 305655
#   ymin: 1718895
#
# conus-wkt: 'PROJCS["Albers",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378140,298.2569999999957,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4326"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
