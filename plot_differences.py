import geopandas as gpd
import pandas as pd
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/usuario/OneDrive/geo_data")
from geo_funciones import crop_to_bounds


scenarios= ["escenario_1","escenario_2","escenario_3","escenario_4","escenario_5","escenario_6","escenario_7","escenario_8","escenario_9","escenario_10"]
lulcs=["built_land","crops","forest","scrub","grassland","rocks","dunes","river_bank","marshes","inrertidal_flats"]

#load the area of interest
aoi=gpd.read_file("concellos_costeiros.shp")
aoi.to_crs("EPSG:3035",inplace=True)
aoi=aoi.dissolve(as_index=False)




df_totals=pd.DataFrame(columns=lulcs,index=scenarios)
print(df_totals)

#baseline scenarios
for lulc in lulcs:
    print(lulc)
    gdf=gpd.read_file("./data/LULC/"+lulc+".shp")
    gdf=gpd.clip(gdf,aoi)
    df_totals.loc["baseline",lulc]=gdf.area.sum()/(10**6)


for scenario in scenarios:
    print(scenario)
    for lulc in lulcs:
        print(lulc)
        gdf=gpd.read_file("./data/inundaciones/"+scenario+"/LULC/"+lulc+".shp")
        gdf=gpd.clip(gdf,aoi)
        df_totals.loc[scenario,lulc]=gdf.area.sum()/(10**6)

print(df_totals)
df_totals.to_csv("inundaciones.csv")