import geopandas as gpd
import time
t0=time.time()
escenario="escenario_1"
gdf=gpd.read_file("data/inundaciones/"+escenario+"/coastline"+escenario+".shp")
gdf=gdf[gdf["raster_val"]==0]
gdf=gdf.dissolve(by="raster_val",as_index=False)
gdf.to_file("./data/inundaciones/escenario_1/cambiar_nombre.shp")
print(time.time()-t0)