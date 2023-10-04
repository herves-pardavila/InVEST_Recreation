import pandas as pd
import geopandas as gpd
import numpy as np
import sys
import rasterio
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from tqdm import tqdm
sys.path.append("/home/usuario/OneDrive/geo_data")
from geo_funciones import rasterize_vector_with_template
from geo_funciones import create_square_grid,polygonize
from qgis_commands import voronoi_polygons,rasterize,reclass
#from esgf_data import ESGF
from mpl_toolkits.mplot3d import axes3d


def station_means(input_file,output_file):
    df=pd.read_csv(input_file)
    df["Value"].replace(-9999.0,np.nan,inplace=True)
    df["Value"].replace(-9999,np.nan,inplace=True)
    df.loc[df["CodeValidation"].isin([0,2,3,9]),"Value"]=np.nan
    newdf=df.groupby(by="idStation",as_index=False).mean(numeric_only=True)
    gdf=gpd.GeoDataFrame(newdf,crs="EPSG:4326",geometry=gpd.points_from_xy(newdf.lon,newdf.lat))
    gdf.to_crs("EPSG:3035",inplace=True)
    gdf[["idStation","Value","altitude","geometry"]].to_file(output_file)

    return
def set_vector_projection(file,crs):
    
    gdf=gpd.read_file(file)
    gdf.to_crs(crs,inplace=True)
    gdf.to_file(file,index=False)
    return


if __name__ == "__main__":

    #area of interest
    create_square_grid(input_file="./data/dissolved_aoi_only_land.shp",size=1000,output_file="./data/aoi_1km_grid.shp")

    #area of the gridded cells
    aoi=gpd.read_file("./data/aoi_1km_grid.shp")
    aoi["area (km2)"]=aoi.geometry.area/1e6
    aoi.to_file("./data/aoi_1km_grid.shp",index=False)

    #raster for the areaof the aoi cells
    rasterize(extent="-10.0,-6.0,41.5,44 [EPSG:3035]",field="area(km2)",height=1000,width=1000,input_vector="./data/aoi_1km_grid.shp",output_raster="./data/area_cells_1kmgrid.tif",units=0,data_type=6)

    #means of the environmental data
    #station_means(input_file="./data/meteogalicia_data/monthlyNDPP1_RECUENTO_1.5m2010-2019.csv",output_file="./data/meteogalicia_data/meanNDPP1_RECUENTO_1.5m.shp")
    station_means(input_file="./data/meteogalicia_data/dailyHSOL_SUM_1.5m2010-2019.csv",output_file="./data/meteogalicia_data/meanHSOL_SUM_1.5m.shp")
    station_means(input_file="./data/meteogalicia_data/dailyPRED_AVG_1.5m2010-2019.csv",output_file="./data/meteogalicia_data/meanPRED_AVG_1.5m.shp")
    station_means(input_file="./data/meteogalicia_data/dailyTA_AVG_1.5m2010-2019.csv",output_file="./data/meteogalicia_data/meanTA_AVG_1.5m.shp")
    station_means(input_file="./data/meteogalicia_data/dailyVV_AVG_2m2010-2019.csv",output_file="./data/meteogalicia_data/meanVV_AVG_2m.shp")

    #Voronoy polygons
    #voronoi_polygons(5,"./data/meteogalicia_data/meanNDPP1_RECUENTO_1.5m.shp","/./data/meteogalicia_data/voronoiNDPP1_RECUENTO_1.5m.shp")
    voronoi_polygons(5,"./data/meteogalicia_data/meanHSOL_SUM_1.5m.shp","./data/meteogalicia_data/voronoiHSOL_SUM_1.5m.shp")
    voronoi_polygons(5,"./data/meteogalicia_data/meanPRED_AVG_1.5m.shp","./data/meteogalicia_data/voronoiPRED_AVG_1.5m.shp")
    voronoi_polygons(5,"./data/meteogalicia_data/meanTA_AVG_1.5m.shp","./data/meteogalicia_data/voronoiTA_AVG_1.5m.shp")
    voronoi_polygons(15,"./data/meteogalicia_data/meanVV_AVG_2m.shp","./data/meteogalicia_data/voronoiVV_AVG_2m.shp")
    #voronoi_polygons(5,"./data/meteogalicia_data/sst.shp","./data/meteogalicia_data/voronoisst.shp")

    #create raster
    #rasterize_vector_with_template("./data/meteogalicia_data/voronoiPP_SUM_1.5m.shp","./data/raster_template.tif","Value","EPSG:3035",np.float,"./data/meteogalicia_data/voronoiPP_SUM_1.5m.tif")
    rasterize_vector_with_template("./data/meteogalicia_data/voronoiHSOL_SUM_1.5m.shp","./data/raster_template.tif","Value","EPSG:3035",np.float64,"./data/meteogalicia_data/voronoiHSOL_SUM_1.5m.tif")
    rasterize_vector_with_template("./data/meteogalicia_data/voronoiPRED_AVG_1.5m.shp","./data/raster_template.tif","Value","EPSG:3035",np.float64,"./data/meteogalicia_data/voronoiPRED_AVG_1.5m.tif")
    rasterize_vector_with_template("./data/meteogalicia_data/voronoiTA_AVG_1.5m.shp","./data/raster_template.tif","Value","EPSG:3035",np.float64,"./data/meteogalicia_data/voronoiTA_AVG_1.5m.tif")
    rasterize_vector_with_template("./data/meteogalicia_data/voronoiVV_AVG_2m.shp","./data/raster_template.tif","Value","EPSG:3035",np.float64,"./data/meteogalicia_data/voronoIVV_AVG_2m.tif")
    #rasterize_vector_with_template("./data/meteogalicia_data/voronoisst.shp","/home/usuario/OneDrive/recreation/InVEST/raster_template.tif","value","EPSG:3035","/home/usuario/OneDrive/recreation/InVEST/voronoisst.tif")
    
    #Select protected zones only
    gdf=gpd.read_file("./data/protected_zones.shp")
    newgdf=gdf[gdf.protected_==255]
    newgdf.to_file("./data/protected_zones_only.shp",index=False)

    #set vector projection of some layers
    set_vector_projection("./data/restaurants_galicia.shp","EPSG:3035")
    set_vector_projection("./data/cities_galicia.shp","EPSG:3035")
    set_vector_projection("./data/towns_galicia.shp","EPSG:3035")
    set_vector_projection("./data/viewpoints_galicia.shp","EPSG:3035")
    set_vector_projection("./data/lighthouse_galicia.shp","EPSG:3035")
    set_vector_projection("./data/OSMtourism.shp","EPSG:3035")
    set_vector_projection("./data/beaches_galicia.shp","EPSG:3035")
    set_vector_projection("./data/nautic_ports_OSM.shp","EPSG:3035")
    set_vector_projection("./data/trial_OSM.shp","EPSG:3035")
    set_vector_projection("./data/bicycle_routesOSM.shp","EPSG:3035")
    set_vector_projection("./data/railwayOSM.shp","EPSG:3035")

    #rasterize beaches
    gdf=gpd.read_file("./data/playas_GEAMA.shp")
    gdf["Area (m2)"]=gdf.geometry.area
    gdf.to_file("./data/beach_GEAMA.shp",index=False)
    rasterize_vector_with_template("./data/beach_GEAMA.shp","./data/raster_template.tif","Area (m2)","EPSG:3035",np.float64,"./data/beach_GEAMA.tif")

    # #rasterize LULC types 2012
    # rasterize_vector_with_template("./data/LULC/CZ_2018_galicia.shp","./data/raster_template.tif","CODE_5_18","EPSG:3035",np.float64,"./data/LULC/CZ_2018_galicia.tif")
    # #reclassify these lulc 
    # reclass(input_file="./data/LULC/CZ_2018_galicia.tif",rules_file="./reclass_rules.txt",output_file="./data/LULC/CZ_2018_galicia_reclass.tif")

    #vectorize reclassified lulc map ----> I did it in QGIS

    #make separate vector files for each lulc type

    gdf=gpd.read_file("./data/LULC/CZ_2012_galicia_reclass.shp")
    print(gdf)
    gdf[gdf.LULC==1].to_file("./data/LULC/built_land.shp",index=False)
    gdf[gdf.LULC==2].to_file("./data/LULC/arable_land.shp",index=False)
    gdf[gdf.LULC==3].to_file("./data/LULC/crops.shp",index=False)
    gdf[gdf.LULC==4].to_file("./data/LULC/forest.shp",index=False)
    gdf[gdf.LULC==5].to_file("./data/LULC/scrub.shp",index=False)
    gdf[gdf.LULC==6].to_file("./data/LULC/grassland.shp",index=False)
    gdf[gdf.LULC==7].to_file("./data/LULC/rocks.shp",index=False)
    gdf[gdf.LULC==8].to_file("./data/LULC/dunes.shp",index=False)
    gdf[gdf.LULC==9].to_file("./data/LULC/river_bank.shp",index=False)
    gdf[gdf.LULC==10].to_file("./data/LULC/marshes.shp",index=False)
    gdf[gdf.LULC==11].to_file("./data/LULC/inrertidal_flats.shp",index=False)
    gdf[gdf.LULC==12].to_file("./data/LULC/water.shp",index=False)
    gdf[gdf.LULC==13].to_file("./data/LULC/reservoirs.shp",index=False)
    gdf[gdf.LULC==14].to_file("./data/LULC/damaged_area.shp",index=False)

    # ===================================SCENARIO FILES==================================================

    
    #Surface loss of beaches due to sea level rise

    scenarios= ["escenario_1","escenario_2","escenario_3","escenario_4","escenario_5","escenario_6","escenario_7","escenario_8","escenario_9","escenario_10"]
    
    for scenario in scenarios:
        print(scenario)
        floods=gpd.read_file("./data/inundaciones/"+scenario+"/coastline"+scenario+".shp")
        floods.to_crs("EPSG:3035",inplace=True)
        floods=floods[floods.raster_val==1]
        
        # #flooding of the lulc types
        # lulcs=["built_land","crops","forest","scrub","grassland","rocks","dunes","river_bank","marshes","inrertidal_flats"]
        # for lulc in lulcs:
        #     print(lulc)
        #     lista=[]
        #     gdf=gpd.read_file("./data/LULC/%s.shp" %lulc)
        #     gdf=gpd.GeoDataFrame({"fid":[lulc],"geometry":[gdf.unary_union]},crs=gdf.crs)
        #     for i in tqdm(range(1,len(floods))):
        #         lista+=[gdf.intersection(floods.iloc[i].geometry)[0]]
            
        #     newgdf=gpd.GeoDataFrame({"LULC":[lulc]*len(lista)},crs=gdf.crs,geometry=lista)
        #     newgdf=newgdf[newgdf.geometry != Polygon([])]
        #     newgdf.to_file("./data/inundaciones/%s/LULC/%s.shp"%(scenario,lulc),index=False)
        #protected zones    
        gdf=gpd.read_file("./data/protected_zones_only.shp",index=False)
        gdf=gpd.GeoDataFrame({"fid":["Protected zones"],"geometry":[gdf.unary_union]},crs=gdf.crs)
        lista=[]
        for i in tqdm(range(1,len(floods))):
            if gdf.intersection(floods.iloc[i].geometry)[0]!=Polygon([]):
                lista+=[gdf.intersection(floods.iloc[i].geometry)[0]]
        newgdf=gpd.GeoDataFrame({"fid":["Protected zone"]*len(lista)},crs=gdf.crs,geometry=lista)
        newgdf.to_file("./data/inundaciones/%s/protected_zones_only.shp"%(scenario),index=False)

        #flooding of the nautic ports
        gdf=gpd.read_file("./data/nautic_ports_OSM.shp",index=False)
        gdf=gpd.GeoDataFrame({"fid":["Nautic Ports"],"geometry":[gdf.unary_union]},crs=gdf.crs)
        for i in tqdm(range(1,len(floods))):
            if gdf.intersection(floods.iloc[i].geometry)[0]!=Polygon([]):
                lista+=[gdf.intersection(floods.iloc[i].geometry)[0]]
        newgdf=gpd.GeoDataFrame({"fid":["Nautic Port"]*len(lista)},crs=gdf.crs,geometry=lista)
        newgdf.to_file("./data/inundaciones/%s/nautic_ports_OSM.shp"%(scenario),index=False)

        #flooding of  tourist infrastructure
        gdf=gpd.read_file("./data/OSMtourism.shp")
        overlay=gpd.overlay(gdf,floods)
        overlay.to_file("./data/inundaciones/%s/OSMtourism.shp"%(scenario),index=False)

        #flooding of  tourist infrastructure
        gdf=gpd.read_file("./data/restaurants_galicia.shp")
        overlay=gpd.overlay(gdf,floods)
        overlay.to_file("./data/inundaciones/%s/restaurants_galicia.shp"%(scenario),index=False)







    #================================== REPRESENTATION 3D===============================================
    # x=np.arange(len(df.LULC.unique()))
    # y=np.arange(len(df.Scenario.unique()))
    # xx,yy=np.meshgrid(x,y)
    # x,y=xx.ravel(),yy.ravel()

    # fig = plt.figure(figsize=(8, 3))
    # ax = fig.add_subplot(111, projection='3d')
    # fig.subplots_adjust(bottom=0.2)
    # colors=["blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","red","red","red","red","red","red","red","red","red","red","green","green","green","green","green","green","green","green","green","green","pink","pink","pink","pink","pink","pink","pink","pink","pink","pink"]
    # ax.bar3d(x,y,np.zeros_like(df["Initial_area (km2)"]),0.5,0.1,df["Initial_area (km2)"],color=colors)
    # # for xs, ys, zs, label in zip(x,y,df["Area_loss (km2)"],df["Area_loss (%)"]):
    # #     my_label=str(label)[0:4]
    # #     ax.text(xs,ys,zs,my_label,fontsize=18)
    # ax.set_yticklabels([""])
    # ax.text(11,0,-1,scenarios[0],fontsize=14)
    # ax.text(11,1,-1,scenarios[1],fontsize=14)
    # ax.text(11,2,-1,scenarios[2],fontsize=14)
    # ax.text(11,3,-1,scenarios[-1],fontsize=14)
    # ax.set_xticks(x)
    # ax.set_yticks(y)    
    # ax.tick_params(axis='z', which='major', labelsize=14)
    # ax.set_xticklabels(["Construcción", "Bosque", "Matorral","Rocas","Marismas","Hierba", "Intermareal","Cultivable","Cultivo","Puerto"])
    # ax.set_zlabel(r"Superficie (km$^2$)",fontsize=14)
    # #plt.zticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # plt.xticks(fontsize=14,rotation = 45)

    # plt.show()

    #================================== REPRESENTATION 3D===============================================

    # fig1=plt.figure()
    # ax1=fig1.add_subplot(111)
    # fig1.subplots_adjust(bottom=0.3)
    # ax1.set_ylabel(r"Superficie (km$^2$)")
    # ax1.set_xticklabels(["Construcción", "Bosque", "Matorral","Rocas","Marismas","Hierba", "Intermareal","Cultivable","Cultivo","Puerto"])
    # ax1.bar(lulcs,df[df.Scenario == scenarios[0]]["Initial_area (km2)"])
    # plt.xticks(fontsize=14,rotation = 45)
    # plt.show()
   
    
    # fig2=plt.figure()
    # ax2=fig2.add_subplot(111)
    # ax2.set_ylabel("Loss relative to their area without climate change")
    # ax2.bar(lulcs,np.array(losses)/np.array(areas))
    # fig2.savefig("/home/usuario/OneDrive/recreation/InVEST/lulc/relative_losses_%s.png"%scenario)
    # plt.show()
    # ==========================EVIRONMNETAL CLIMATE CHANGE PREDICTIONS====================
    #con=ESGF("/home/usuario/OneDrive/geo_data/CMIP6/HighResMIP/ps_Amon_HiRAM-SIT-HR_hist-1950_r1i1p1f1_gn_201401-201412.nc")