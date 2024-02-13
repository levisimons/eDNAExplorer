import ee

class VectorDataset:
    
    def __init__(self, snippet, property: str, name=None, map_params={}):
        """
        Represent a Google Earth Engine vector dataset and the information we need to use it.
        
        snippet: The ee.FeatureCollection object representing the dataset.
            This is the 'Earth Engine snippet' displayed on the dataset's page in the Earth 
            Engine catalog. e.g.: ee.FeatureCollection('EPA/Ecoregions/2013/L3')
        property: The key of the desired feature property
        name: The column name to display for the data gathered from this dataset.
        map_params: Settings for how to display the data e.g. palette, min, max
        """
        self.data = snippet
        self.property = property
        self.name = name or property
        self.map_params = map_params
        
    # Note: This function is an argument to map(). Arguments to map() cannot print anything
    # or call getInfo(). Doing so results in an EEException: ValueNode empty
    # source: https://gis.stackexchange.com/questions/345598/mapping-simp>le-function-to-print-date-and-time-stamp-on-google-earth-engine-pyth
    def get_sample_area_data(self, sample_area: ee.Feature) -> str:
        """
        Return the value from the dataset to assign to the sample area.
        """
        # Get a FeatureCollection storing the overlaps between the sample area and the dataset
        overlaps = self.data.filterBounds(sample_area.geometry()).map(
            lambda feature: feature.intersection(sample_area.geometry())
        )

        result = ee.Algorithms.If(
            # If there is exactly 1 overlapping dataset feature, return its value
            overlaps.size().eq(1),
            sample_area.set(self.name, overlaps.first().get(self.property)),

            ee.Algorithms.If(
                # If there are 0 overlapping dataset features, return the null value
                overlaps.size().eq(0),
                sample_area.set({}),
                #sample_area.set(self.name, 
                #                self.get_nearest_feature(sample_area).get(self.property)),

                # Otherwise, there must be >1 features overlapping the sample area
                # Return the value of the one with the largest overlap
                sample_area.set(self.name, 
                                self.get_predominant_feature(overlaps).get(self.property))
            )
        )


        return result

    def get_nearest_feature(self, sample_area: ee.Feature) -> str:
        """
        To be used when the sample area doesn't overlap the dataset at all.
        
        Get the dataset feature that is nearest to the sample area, and
        return the value of its dataset.property.
        
        param sample_area: ee.Feature representing the sample area of interest
        """

        # Define a filter to get all dataset features within 10000 meters of the sample area
        spatialFilter = ee.Filter.withinDistance(
            distance=10000,
            leftField='.geo',
            rightField='.geo',
            maxError=10
        )
        # Define a join that will return only the 'best' (nearest) match
        saveBestJoin = ee.Join.saveBest(
          matchKey='closestFeature',
          measureKey='distance'
        )
        # Apply the join, using the distance filter to define match quality
        # Get the only feature in the resulting FeatureCollection
        result = ee.Feature(saveBestJoin.apply(
            ee.FeatureCollection(sample_area),
            self.data,
            spatialFilter
        ).first())

        # Return the closest dataset feature
        return ee.Feature(result.get('closestFeature'))
    
    def get_predominant_feature(self, overlaps: ee.FeatureCollection) -> str:
        """
        To be used when the sample area overlaps more than one dataset feature.
        
        Return the value of 'property' for the largest overlap.
        """
        # Add 'area' as a property to each feature. This is the area in square meters 
        # of the intersection of the ecoregion feature and the sample area.
        overlaps = overlaps.map(
            lambda feature: feature.set({'area': feature.geometry().area()}))

        # Find the maximum area among all the overlaps
        max_area = overlaps.aggregate_max('area')

        # Return the overlap with the largest area
        return ee.Feature(overlaps.filter(ee.Filter.gte('area', max_area)).first())