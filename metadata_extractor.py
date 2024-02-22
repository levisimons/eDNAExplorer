import ee
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
from ee.ee_exception import EEException
from RasterDataset import RasterDataset
from VectorDataset import VectorDataset
import csv
import argparse
import sys

# script arguments
parser = argparse.ArgumentParser(description='Google Earth Engine metadata extractor')
parser.add_argument('--input', type=str, help='Header row of InputMetadata.csv')
parser.add_argument('--output', type=str, help='Path to output file')
parser.add_argument('--account', type=str, help='Service Account name')
parser.add_argument('--credentials', type=str, help='Path to Service Account credentials')
args = parser.parse_args()

outfile = args.output
csv_file = args.input # split by column
service_account = args.account
service_credentials = args.credentials

def mask_clouds(image: ee.Image) -> ee.Image:
  #Return the image with clouds masked
  qa = image.select('QA60')
  
  # Bits 10 and 11 are clouds and cirrus, respectively.
  cloudBitMask = int(2 ** 10)
  cirrusBitMask = int(2 ** 11)
  
  # Both flags should be set to zero, indicating clear conditions.
  mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))
  image = image.updateMask(mask)
  
  return image


def filter_StudyArea(collection: ee.ImageCollection) -> ee.ImageCollection:
  #Filter the collection to just those images overlapping the study area surrounding the sample coordinates'''
  #Define a rectangular area surrounding the sample coordinates with a +/- 0.5 degree buffer.
  StudyArea = ee.Geometry.Polygon([
      ee.Geometry.Point(WestBound, SouthBound),
      ee.Geometry.Point(WestBound, NorthBound),
      ee.Geometry.Point(EastBound, NorthBound),
      ee.Geometry.Point(EastBound, SouthBound)
    ])
  # Return a collection of just those images that are within the study area
  return collection.filterBounds(StudyArea)

# Here is an extra function for reprojecting an Image to a different scale. 
# Might be useful for something, not sure what effect this has exactly

def reproject(image: ee.Image, scale: int) -> ee.Image:
  #Resample and reproject an image to a different resolution.
  #param image: Image to reproject
  #param scale: Desired pixel width in meters
  return image.resample('bilinear').reproject(  # Use bilinear interpolation method
    crs=image.projection().crs(),  # Keep the same map projection
    scale=scale  # Change the scale
  )

# Set a minimum value for Spatial_Uncertainty
def set_min_spatial_uncertainty(value):
    if pd.isna(value) or value < 30:
        return 30.0
    return value

# Read the relevant columns from the CSV into a dataframe
samples = pd.read_csv(csv_file, usecols=['Sample ID', 'Longitude', 'Latitude', 'Sample Date', 'Spatial Uncertainty'])
samples.rename(columns={'Sample ID': 'Sample_ID', 'Sample Date': 'Sample_Date', 'Spatial Uncertainty': 'Spatial_Uncertainty'}, inplace=True)
samples = samples.astype({'Sample_ID': 'object', 'Longitude': 'float64', 'Latitude': 'float64', 'Sample_Date': 'object', 'Spatial_Uncertainty': 'float64'})

# Apply set_min_spatial_uncertainty to 'Spatial_Uncertainty' column
samples['Spatial_Uncertainty'] = samples['Spatial_Uncertainty'].apply(set_min_spatial_uncertainty)

# Calculate a buffer from the Spatial_Uncertainty in meters.
BufferRadius = samples['Spatial_Uncertainty'].values[0]

# Convert missing spatial uncertainties to a default radius in meters.
# Reset buffer if it is too small.
if pd.isna(BufferRadius):
    BufferRadius = 30.0
if BufferRadius < 30:
    BufferRadius = 30.0

# Convert missing spatial uncertainties to a default radius in meters.
#samples['Spatial_Uncertainty'] = samples['Spatial_Uncertainty'].fillna(30)
#if not str(BufferRadius).isnumeric():
#  BufferRadius = 30
#  samples['Spatial_Uncertainty'].values[0] = BufferRadius

# Reset buffer if it is too small.
#samples['Spatial_Uncertainty'] = samples['Spatial_Uncertainty'].apply(lambda x: max(x, 30))
#if BufferRadius < 3:
#  BufferRadius = 3

# Remove entries with missing data
samples = samples.dropna()
if(samples.empty):
    print("Samples empty. Exiting")
    sys.exit()

credentials = ee.ServiceAccountCredentials(service_account, service_credentials)
ee.Initialize(credentials)

#Read in Sample_Date as a date.
samples['Sample_Date'] = pd.to_datetime(samples['Sample_Date'],infer_datetime_format=True)
#Define earliest and latest times using Sample_Dates with a 6 month leading buffer.
EarliestTime = min(samples['Sample_Date']) - relativedelta(months=6)
LatestTime = max(samples['Sample_Date'])

#Define geographic bounds using sample coordinates and a +/- 0.5 degree buffer.
EastBound = max(samples['Longitude']) + 0.5
WestBound = min(samples['Longitude']) - 0.5
SouthBound = min(samples['Latitude']) - 0.5
NorthBound = max(samples['Latitude']) + 0.5

# Define a Feature for each sample area within the Spatial_Uncertainty of each sample location.
sample_areas = []
for sample in samples.itertuples():
  # Store the important data as properties of the feature
  sample_areas.append(
      ee.Feature(ee.Geometry.Point(sample.Longitude, sample.Latitude).buffer(BufferRadius))
      .set('name', sample.Sample_ID)
      .set('Sample_Date', sample.Sample_Date)
      .set('Longitude', sample.Longitude)
      .set('Latitude', sample.Latitude)
      .set('Spatial_Uncertainty', sample.Spatial_Uncertainty))

# Define a FeatureCollection containing all the sample areas
sample_areas = ee.FeatureCollection(sample_areas)

#Get bioclim data.
bioclim = [
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)),
        band='bio01'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)),
        band='bio02'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'), 
        band='bio03'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(100)), 
        band='bio04'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)), 
        band='bio05'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)), 
        band='bio06'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)), 
        band='bio07'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)), 
        band='bio08'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)), 
        band='bio09'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)), 
        band='bio10'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO').divide(ee.Image(10)), 
        band='bio11'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'),
        band='bio12'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'), 
        band='bio13'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'), 
        band='bio14'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'),
        band='bio15'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'),
        band='bio16'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'), 
        band='bio17'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'), 
        band='bio18'
    ),
    RasterDataset(
        snippet=ee.Image('WORLDCLIM/V1/BIO'), 
        band='bio19'
    )
]


soil = [
    RasterDataset(  # Soil taxonomy great groups, a classification number from 0 to 433
        snippet=ee.Image('OpenLandMap/SOL/SOL_GRTGROUP_USDA-SOILTAX_C/v01'),
        band='grtgroup',
        categorical=True
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').divide(ee.Image(10)),
        band='b0',
        name='soil pH in H20 at 0 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').divide(ee.Image(10)),
        band='b10',
        name='soil pH in H20 at 10 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').divide(ee.Image(10)),
        band='b30',
        name='soil pH in H20 at 30 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').divide(ee.Image(10)),
        band='b60',
        name='soil pH in H20 at 60 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').divide(ee.Image(10)),
        band='b100',
        name='soil pH in H20 at 100 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02').divide(ee.Image(10)),
        band='b200',
        name='soil pH in H20 at 200 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').divide(ee.Image(5)),
        band='b0',
        name='soil organic carbon at 0 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').divide(ee.Image(5)),
        band='b10',
        name='soil organic carbon at 10 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').divide(ee.Image(5)),
        band='b30',
        name='soil organic carbon at 30 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').divide(ee.Image(5)),
        band='b60',
        name='soil organic carbon at 60 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').divide(ee.Image(5)),
        band='b100',
        name='soil organic carbon at 100 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').divide(ee.Image(5)),
        band='b200',
        name='soil organic carbon at 200 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b0',
        name='soil sand at 0 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b10',
        name='soil sand at 10 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b30',
        name='soil sand at 30 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b60',
        name='soil sand at 60 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b100',
        name='soil sand at 100 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b200',
        name='soil sand at 200 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b0',
        name='soil clay at 0 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b10',
        name='soil clay at 10 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b30',
        name='soil clay at 30 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b60',
        name='soil clay at 60 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b100',
        name='soil clay at 100 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02'),
        band='b200',
        name='soil clay at 200 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02').divide(ee.Image(10)),
        band='b0',
        name='soil bulk density at 0 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02').divide(ee.Image(10)),
        band='b10',
        name='soil bulk density at 10 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02').divide(ee.Image(10)),
        band='b30',
        name='soil bulk density at 30 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02').divide(ee.Image(10)),
        band='b60',
        name='soil bulk density at 60 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02').divide(ee.Image(10)),
        band='b100',
        name='soil bulk density at 100 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_BULKDENS-FINEEARTH_USDA-4A1H_M/v02').divide(ee.Image(10)),
        band='b200',
        name='soil bulk density at 200 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02'),
        band='b0',
        name='soil texture class at 0 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02'),
        band='b10',
        name='soil texture class at 10 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02'),
        band='b30',
        name='soil texture class at 30 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02'),
        band='b60',
        name='soil texture class at 60 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02'),
        band='b100',
        name='soil texture class at 100 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02'),
        band='b200',
        name='soil texture class at 200 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01'),
        band='b0',
        name='Percent volume soil water content at 0 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01'),
        band='b10',
        name='Percent volume soil water content at 10 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01'),
        band='b30',
        name='Percent volume soil water content at 30 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01'),
        band='b60',
        name='Percent volume soil water content at 60 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01'),
        band='b100',
        name='Percent volume soil water content at 100 cm'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01'),
        band='b200',
        name='Percent volume soil water content at 200 cm'
    )
]


terrain = [
    RasterDataset(
        snippet=ee.Image('CGIAR/SRTM90_V4'),
        band='elevation'
    ),
    RasterDataset(
        snippet=ee.Terrain.slope(ee.Image('CGIAR/SRTM90_V4')),
        band='slope'
    ),
    RasterDataset(
        snippet=ee.Terrain.aspect(ee.Image('CGIAR/SRTM90_V4')),
        band='aspect'
    ),
]


HII = [
    RasterDataset(
        snippet=ee.ImageCollection('CSP/HM/GlobalHumanModification'),
        band='gHM'
    ),
]


landsat = [
    RasterDataset(
        snippet=ee.ImageCollection('LANDSAT/LC08/C01/T1_32DAY_NDVI').filterDate(
            EarliestTime, LatestTime
        ),
        band='NDVI',
        map_params={'min': 0, 'max': 1, 'palette': [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
    '012E01', '011D01', '011301'
  ]}
    ),
    RasterDataset(
        snippet=ee.ImageCollection('LANDSAT/LC08/C01/T1_8DAY_EVI').filterDate(
            EarliestTime, LatestTime
        ),
        band='EVI'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('LANDSAT/LC08/C01/T1_8DAY_NBRT').filterDate(
            EarliestTime, LatestTime
        ),
        band='NBRT'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('LANDSAT/LC08/C01/T1_ANNUAL_GREENEST_TOA').filterDate(
            EarliestTime, LatestTime
        ),
        band='greenness'
    )
]


sentinel_preprocessed = filter_StudyArea(
                            ee.ImageCollection("COPERNICUS/S2")
                        ).filterDate(
                            EarliestTime,
                            LatestTime
                        ).filter(
                            ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)
                        ).map(
                            lambda img: mask_clouds(img).divide(
                                ee.Image(10000)
                            )
                        ).median()

sentinel = [
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B1',
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B2',
        map_params={'min': 0, 'max': 1}
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B3',
        map_params={'min': 0, 'max': 1}
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B4',
        map_params={'min': 0, 'max': 1}
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B5'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B6'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B7'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B8'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B8A'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B9'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B10'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B11'
    ),
    RasterDataset(
        snippet=sentinel_preprocessed,
        band='B12'
    ),
]


Population = [
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='population',
        name='total population per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_0',
        name='total male population age 0 to 1 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_1',
        name='total male population age 1 to 4 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_5',
        name='total male population age 5 to 9 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_10',
        name='total male population age 10 to 14 per hectare'
    ),
     RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_15',
        name='total male population age 15 to 19 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_20',
        name='total male population age 20 to 24 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_25',
        name='total male population age 25 to 29 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_30',
        name='total male population age 30 to 34 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_35',
        name='total male population age 35 to 39 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_40',
        name='total male population age 40 to 44 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_45',
        name='total male population age 45 to 49 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_50',
        name='total male population age 50 to 54 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_55',
        name='total male population age 55 to 59 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_60',
        name='total male population age 60 to 64 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_65',
        name='total male population age 65 to 65 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_70',
        name='total male population age 70 to 74 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_75',
        name='total male population age 75 to 79 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='M_80',
        name='total male population age 80 to 84 per hectare'
    ),
        RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_0',
        name='total female population age 0 to 1 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_1',
        name='total female population age 1 to 4 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_5',
        name='total female population age 5 to 9 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_10',
        name='total female population age 10 to 14 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_15',
        name='total female population age 15 to 19 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_20',
        name='total female population age 20 to 24 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_25',
        name='total female population age 25 to 29 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_30',
        name='total female population age 30 to 34 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_35',
        name='total female population age 35 to 39 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_40',
        name='total female population age 40 to 44 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_45',
        name='total female population age 45 to 49 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_50',
        name='total female population age 50 to 54 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_55',
        name='total female population age 55 to 59 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_60',
        name='total female population age 60 to 64 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_65',
        name='total female population age 65 to 65 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_70',
        name='total female population age 70 to 74 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_75',
        name='total female population age 75 to 79 per hectare'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('WorldPop/GP/100m/pop_age_sex_cons_unadj'),
        band='F_80',
        name='total female population age 80 to 84 per hectare'
    )
]



PotentialBiomes = [
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_BIOME-TYPE_BIOME00K_C/v01'),
        band='biome_type'
    )
]



HealthCareAccess = [
    RasterDataset(
        snippet=ee.Image('Oxford/MAP/accessibility_to_healthcare_2019'),
        band='accessibility',
        name='Travel time to nearest healthcare facility'
    ),
    RasterDataset(
        snippet=ee.Image('Oxford/MAP/accessibility_to_healthcare_2019'),
        band='accessibility_walking_only',
        name='Waling time to nearest healthcare facility'
    )
]



GlobalFriction = [
    RasterDataset(
        snippet=ee.Image('Oxford/MAP/friction_surface_2019'),
        band='friction',
        name='Terrestrial travel speed'
    ),
    RasterDataset(
        snippet=ee.Image('Oxford/MAP/friction_surface_2019'),
        band='friction_walking_only',
        name='Walking speed'
    )
]



LightPollution = [
    RasterDataset(
        snippet=ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG').filterDate(
            EarliestTime, LatestTime
        ),
        band='avg_rad',
        name='Average_Radiance'
    )
]



CHILI = [
    RasterDataset(
        snippet=ee.Image('CSP/ERGo/1_0/Global/SRTM_CHILI'),
        band='constant',
        name='Continuous Heat Insolation Load Index'
    )
]



Landform = [
    RasterDataset(
        snippet=ee.Image('CSP/ERGo/1_0/Global/SRTM_landforms'),
        band='constant',
        name='Landform'
    )
]



FAPAR = [
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='jan',
        name='FAPAR January'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='feb',
        name='FAPAR February'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='mar',
        name='FAPAR March'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='apr',
        name='FAPAR April'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='may',
        name='FAPAR May'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='jun',
        name='FAPAR June'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='jul',
        name='FAPAR July'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='aug',
        name='FAPAR August'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='sep',
        name='FAPAR September'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='oct',
        name='FAPAR October'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='nov',
        name='FAPAR November'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='dec',
        name='FAPAR December'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='annual',
        name='FAPAR Annual'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/PNV/PNV_FAPAR_PROBA-V_D/v01'),
        band='annualdiff',
        name='FAPAR Annual difference'
    )
]



MonthlyPrecipitation = [
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='jan',
        name='January rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='feb',
        name='February rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='mar',
        name='March rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='apr',
        name='April rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='may',
        name='May rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='jun',
        name='June rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='jul',
        name='July rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='aug',
        name='August rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='sep',
        name='September rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='oct',
        name='October rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='nov',
        name='November rainfall'
    ),
    RasterDataset(
        snippet=ee.Image('OpenLandMap/CLM/CLM_PRECIPITATION_SM2RAIN_M/v01'),
        band='dec',
        name='December rainfall'
    )
]



ecoregions = [
    VectorDataset(
        snippet=ee.FeatureCollection('RESOLVE/ECOREGIONS/2017').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='ECO_NAME'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('RESOLVE/ECOREGIONS/2017').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='REALM'
    )
]



watersheds = [
    VectorDataset(
        snippet=ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_9').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='HYBAS_ID'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_9').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='SUB_AREA'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_9').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='UP_AREA'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_9').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='ENDO'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_9').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='COAST'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_9').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='ORDER'
    )
]



ProtectedAreas = [
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='WDPAID'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='WDPA_PID'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='DESIG_ENG'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='IUCN_CAT'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='MARINE'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='REP_AREA'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='STATUS_YR'
    ),
    VectorDataset(
        snippet=ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(ee.Geometry.Rectangle(WestBound,SouthBound,EastBound,NorthBound)),
        property='GOV_TYPE'
    )
]



OceanColor = [
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='chlor_a',
        name='Chlorophyll a concentration'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='poc',
        name='Particulate organic carbon'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_412',
        name='Remote sensing reflectance at band 412nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_443',
        name='Remote sensing reflectance at band 443nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_469',
        name='Remote sensing reflectance at band 469nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_488',
        name='Remote sensing reflectance at band 488nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_531',
        name='Remote sensing reflectance at band 531nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_547',
        name='Remote sensing reflectance at band 547nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_555',
        name='Remote sensing reflectance at band 555nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_645',
        name='Remote sensing reflectance at band 645nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_667',
        name='Remote sensing reflectance at band 667nm'
    ),
    RasterDataset(
        snippet=ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI').filterDate(
            EarliestTime, LatestTime
        ),
        band='Rrs_678',
        name='Remote sensing reflectance at band 678nm'
    )
]



Fishing = [
    RasterDataset(
        snippet=ee.ImageCollection('GFW/GFF/V1/fishing_hours').filterDate(
            EarliestTime, LatestTime
        ).median(),
        band='drifting_longlines',
    ),
    RasterDataset(
        snippet=ee.ImageCollection('GFW/GFF/V1/fishing_hours').filterDate(
            EarliestTime, LatestTime
        ).median(),
        band='fixed_gear',
    ),
    RasterDataset(
        snippet=ee.ImageCollection('GFW/GFF/V1/fishing_hours').filterDate(
            EarliestTime, LatestTime
        ).median(),
        band='other_fishing',
    ),
    RasterDataset(
        snippet=ee.ImageCollection('GFW/GFF/V1/fishing_hours').filterDate(
            EarliestTime, LatestTime
        ).median(),
        band='purse_seines',
    ),
    RasterDataset(
        snippet=ee.ImageCollection('GFW/GFF/V1/fishing_hours').filterDate(
            EarliestTime, LatestTime
        ).median(),
        band='squid_jigger',
    ),
    RasterDataset(
        snippet=ee.ImageCollection('GFW/GFF/V1/fishing_hours').filterDate(
            EarliestTime, LatestTime
        ).median(),
        band='trawlers',
    )
]



SalinityAndTemperature_preprocessed = ee.ImageCollection('HYCOM/sea_temp_salinity').filterDate(EarliestTime,LatestTime)
#SalinityAndTemperature = SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp__0').multiply(0.001).add(20)
SalinityAndTemperature = [
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_0').multiply(0.001).add(20),
        band='water_temp_0',
        name='Sea water temperature at a depth of 0m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_0').multiply(0.001).add(20),
        band='salinity_0',
        name='salinity units at depth of 0m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_2').multiply(0.001).add(20),
        band='water_temp_2',
        name='Sea water temperature at a depth of 2m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_2').multiply(0.001).add(20),
        band='salinity_2',
        name='salinity units at depth of 2m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_4').multiply(0.001).add(20),
        band='water_temp_4',
        name='Sea water temperature at a depth of 4m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_4').multiply(0.001).add(20),
        band='salinity_4',
        name='salinity units at depth of 4m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_6').multiply(0.001).add(20),
        band='water_temp_6',
        name='Sea water temperature at a depth of 6m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_6').multiply(0.001).add(20),
        band='salinity_6',
        name='salinity units at depth of 6m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_8').multiply(0.001).add(20),
        band='water_temp_8',
        name='Sea water temperature at a depth of 8m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_8').multiply(0.001).add(20),
        band='salinity_8',
        name='salinity units at depth of 8m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_10').multiply(0.001).add(20),
        band='water_temp_10',
        name='Sea water temperature at a depth of 10m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_10').multiply(0.001).add(20),
        band='salinity_10',
        name='salinity units at depth of 10m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_12').multiply(0.001).add(20),
        band='water_temp_12',
        name='Sea water temperature at a depth of 12m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_12').multiply(0.001).add(20),
        band='salinity_12',
        name='salinity units at depth of 12m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_15').multiply(0.001).add(20),
        band='water_temp_15',
        name='Sea water temperature at a depth of 15m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_15').multiply(0.001).add(20),
        band='salinity_15',
        name='salinity units at depth of 15m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_20').multiply(0.001).add(20),
        band='water_temp_20',
        name='Sea water temperature at a depth of 20m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_20').multiply(0.001).add(20),
        band='salinity_20',
        name='salinity units at depth of 20m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_25').multiply(0.001).add(20),
        band='water_temp_25',
        name='Sea water temperature at a depth of 25m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_25').multiply(0.001).add(20),
        band='salinity_25',
        name='salinity units at depth of 25m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_30').multiply(0.001).add(20),
        band='water_temp_30',
        name='Sea water temperature at a depth of 30m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_30').multiply(0.001).add(20),
        band='salinity_30',
        name='salinity units at depth of 30m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_35').multiply(0.001).add(20),
        band='water_temp_35',
        name='Sea water temperature at a depth of 35m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_35').multiply(0.001).add(20),
        band='salinity_35',
        name='salinity units at depth of 35m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_40').multiply(0.001).add(20),
        band='water_temp_40',
        name='Sea water temperature at a depth of 40m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_40').multiply(0.001).add(20),
        band='salinity_40',
        name='salinity units at depth of 40m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_45').multiply(0.001).add(20),
        band='water_temp_45',
        name='Sea water temperature at a depth of 45m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_45').multiply(0.001).add(20),
        band='salinity_45',
        name='salinity units at depth of 45m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_50').multiply(0.001).add(20),
        band='water_temp_50',
        name='Sea water temperature at a depth of 50m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_50').multiply(0.001).add(20),
        band='salinity_50',
        name='salinity units at depth of 50m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_60').multiply(0.001).add(20),
        band='water_temp_60',
        name='Sea water temperature at a depth of 60m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_60').multiply(0.001).add(20),
        band='salinity_60',
        name='salinity units at depth of 60m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_70').multiply(0.001).add(20),
        band='water_temp_70',
        name='Sea water temperature at a depth of 70m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_70').multiply(0.001).add(20),
        band='salinity_70',
        name='salinity units at depth of 70m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_80').multiply(0.001).add(20),
        band='water_temp_80',
        name='Sea water temperature at a depth of 80m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_80').multiply(0.001).add(20),
        band='salinity_80',
        name='salinity units at depth of 80m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_90').multiply(0.001).add(20),
        band='water_temp_90',
        name='Sea water temperature at a depth of 90m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_90').multiply(0.001).add(20),
        band='salinity_90',
        name='salinity units at depth of 90m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_100').multiply(0.001).add(20),
        band='water_temp_100',
        name='Sea water temperature at a depth of 100m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_100').multiply(0.001).add(20),
        band='salinity_100',
        name='salinity units at depth of 100m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_125').multiply(0.001).add(20),
        band='water_temp_125',
        name='Sea water temperature at a depth of 125m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_125').multiply(0.001).add(20),
        band='salinity_125',
        name='salinity units at depth of 125m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_150').multiply(0.001).add(20),
        band='water_temp_150',
        name='Sea water temperature at a depth of 150m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_150').multiply(0.001).add(20),
        band='salinity_150',
        name='salinity units at depth of 150m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_200').multiply(0.001).add(20),
        band='water_temp_200',
        name='Sea water temperature at a depth of 200m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_200').multiply(0.001).add(20),
        band='salinity_200',
        name='salinity units at depth of 200m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_250').multiply(0.001).add(20),
        band='water_temp_250',
        name='Sea water temperature at a depth of 250m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_250').multiply(0.001).add(20),
        band='salinity_250',
        name='salinity units at depth of 250m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_300').multiply(0.001).add(20),
        band='water_temp_300',
        name='Sea water temperature at a depth of 300m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_300').multiply(0.001).add(20),
        band='salinity_300',
        name='salinity units at depth of 300m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_350').multiply(0.001).add(20),
        band='water_temp_350',
        name='Sea water temperature at a depth of 350m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_350').multiply(0.001).add(20),
        band='salinity_350',
        name='salinity units at depth of 350m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_400').multiply(0.001).add(20),
        band='water_temp_400',
        name='Sea water temperature at a depth of 400m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_400').multiply(0.001).add(20),
        band='salinity_400',
        name='salinity units at depth of 400m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_500').multiply(0.001).add(20),
        band='water_temp_500',
        name='Sea water temperature at a depth of 500m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_500').multiply(0.001).add(20),
        band='salinity_500',
        name='salinity units at depth of 500m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_600').multiply(0.001).add(20),
        band='water_temp_600',
        name='Sea water temperature at a depth of 600m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_600').multiply(0.001).add(20),
        band='salinity_600',
        name='salinity units at depth of 600m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_700').multiply(0.001).add(20),
        band='water_temp_700',
        name='Sea water temperature at a depth of 700m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_700').multiply(0.001).add(20),
        band='salinity_700',
        name='salinity units at depth of 700m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_800').multiply(0.001).add(20),
        band='water_temp_800',
        name='Sea water temperature at a depth of 800m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_800').multiply(0.001).add(20),
        band='salinity_800',
        name='salinity units at depth of 800m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_900').multiply(0.001).add(20),
        band='water_temp_900',
        name='Sea water temperature at a depth of 900m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_900').multiply(0.001).add(20),
        band='salinity_900',
        name='salinity units at depth of 900m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_1000').multiply(0.001).add(20),
        band='water_temp_1000',
        name='Sea water temperature at a depth of 1000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_1000').multiply(0.001).add(20),
        band='salinity_1000',
        name='salinity units at depth of 1000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_1250').multiply(0.001).add(20),
        band='water_temp_1250',
        name='Sea water temperature at a depth of 1250m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_1250').multiply(0.001).add(20),
        band='salinity_1250',
        name='salinity units at depth of 1250m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_1500').multiply(0.001).add(20),
        band='water_temp_1500',
        name='Sea water temperature at a depth of 1500m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_1500').multiply(0.001).add(20),
        band='salinity_1500',
        name='salinity units at depth of 1500m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_2000').multiply(0.001).add(20),
        band='water_temp_2000',
        name='Sea water temperature at a depth of 2000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_2000').multiply(0.001).add(20),
        band='salinity_2000',
        name='salinity units at depth of 2000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_2500').multiply(0.001).add(20),
        band='water_temp_2500',
        name='Sea water temperature at a depth of 2500m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_2500').multiply(0.001).add(20),
        band='salinity_2500',
        name='salinity units at depth of 2500m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_3000').multiply(0.001).add(20),
        band='water_temp_3000',
        name='Sea water temperature at a depth of 3000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_3000').multiply(0.001).add(20),
        band='salinity_3000',
        name='salinity units at depth of 3000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_4000').multiply(0.001).add(20),
        band='water_temp_4000',
        name='Sea water temperature at a depth of 4000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_4000').multiply(0.001).add(20),
        band='salinity_4000',
        name='salinity units at depth of 4000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'water_temp_5000').multiply(0.001).add(20),
        band='water_temp_5000',
        name='Sea water temperature at a depth of 5000m'
    ),
    RasterDataset(
        snippet=SalinityAndTemperature_preprocessed.limit(1, 'system:time_start', False).first().select(f'salinity_5000').multiply(0.001).add(20),
        band='salinity_5000',
        name='salinity units at depth of 5000m'
    )
]




# Define the list of datasets from which to retrieve data
datasets = bioclim + soil + terrain + HII + landsat + sentinel + Population + PotentialBiomes + GlobalFriction + LightPollution + CHILI + Landform + FAPAR + MonthlyPrecipitation + ecoregions + watersheds + ProtectedAreas + Fishing + OceanColor + SalinityAndTemperature

for dataset in datasets:
    # print(dataset.name)
    try:
      dataset.data.getInfo()
      sample_areas = sample_areas.map(lambda feature: dataset.get_sample_area_data(feature))
    except EEException:
    #   print(f"dataset {dataset.name} has no input value.")
      sample_areas.map(lambda feature: dataset.get_sample_area_data(feature))

# Get the name of each dataset, which will be the column headers
dataset_names = [dataset.name for dataset in datasets]

# Fetch the table data
table_data = sample_areas.getInfo().get('features', [])

rows = []
#Export extracted row.  
for feature in table_data:
  properties = feature['properties']
  header = ['name', 'Sample_Date', 'Latitude', 'Longitude', 'Spatial_Uncertainty'] + dataset_names
  row = [properties.get(key) for key in header]
  date_value = row[1]['value']
  date = datetime.fromtimestamp(date_value / 1000).strftime("%Y-%m-%d %H:%M:%S")
  row[1] = date
  rows.append(row)

with open(outfile, 'w', newline='') as out:
    writer = csv.writer(out)
    writer.writerow(header)
    for row in rows:
        print(row)
        writer.writerow(row)

print(len(rows))