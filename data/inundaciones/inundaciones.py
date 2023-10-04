import rasterio
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import os
import sys
sys.path.append("/home/usuario/OneDrive/geo_data")
from geo_funciones import polygonize


if __name__== "__main__":

    escenario="escenario_10"
    files=os.listdir(escenario)
    gdf=polygonize("./"+escenario+"/"+files[0])
    crs=gdf.crs
    for file in files[1:]:
        print(file)
        gdfnew=polygonize("./"+escenario+"/"+file)
        gdf=pd.concat([gdf,gdfnew])
    gdf=gpd.GeoDataFrame(data=gdf,crs=crs,geometry=gdf.geometry)
    gdf.to_file("./"+escenario+"/"+"coastline"+escenario+".shp")
    gdf.plot("raster_val",legend=True)
    plt.show()
