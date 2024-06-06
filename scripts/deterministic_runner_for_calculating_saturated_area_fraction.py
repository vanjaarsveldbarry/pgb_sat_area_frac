#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import shutil
import sys
import datetime

import pcraster as pcr

from pcraster.framework import DynamicModel
from pcraster.framework import DynamicFramework

from currTimeStep import ModelTime

import ncConverter_for_discharge_30sec as netcdf_writer
from tools import virtualOS as vos

import logging
logger = logging.getLogger(__name__)


class DeterministicRunner(DynamicModel):

    def __init__(self, modelTime, model_setup):
        DynamicModel.__init__(self)

        # initiate model time
        self.modelTime = modelTime        


        # get the model setup
        self.model_setup = model_setup
        
        # set clone
        self.clone = self.model_setup["clone_file"]
        pcr.setclone(self.clone)
        
        # output and tmp folders
        self.output_folder = model_setup['output_dir']
        self.tmp_folder    = model_setup['tmp_dir']

        # read ldd
        self.ldd = vos.readPCRmapClone(v                = self.model_setup["ldd_file"], \
                                       cloneMapFileName = self.clone, \
                                       tmpDir           = self.tmp_folder, \
                                       absolutePath     = None, \
                                       isLddMap         = True, \
                                       cover            = None, \
                                       isNomMap         = False)
        self.ldd = pcr.lddrepair(pcr.lddrepair(pcr.ldd(self.ldd)))                               

        # read cell area (m2)
        self.cell_area = vos.readPCRmapClone(v                = self.model_setup["cell_area_file"], \
                                             cloneMapFileName = self.clone, \
                                             tmpDir           = self.tmp_folder, \
                                             absolutePath     = None, \
                                             isLddMap         = False, \
                                             cover            = None, \
                                             isNomMap         = False)

        
        # initiate a netcdf writer
        self.netcdf_report = netcdf_writer.PCR2netCDF(self.model_setup["clone_file"])
        self.netcdf_report.createNetCDF(self.model_setup["discharge_output_file"],\
                                        "discharge",\
                                        "m3/s")
        
    def initial(self): 
        
        # read soil parameters from the netCDF file:
        soilParameters = ['resVolWC1',\           
                          'resVolWC2',\
                          'satVolWC1',\
                          'satVolWC2']
        for var in soilParameters:
            soil_input_file = self.model_setup[var]
            vars(self)[var] = vos.readPCRmapClone(soil_input_file, self.cloneMap, self.tmpDir, self.inputDir)
            vars(self)[var] = pcr.scalar(vars(self)[var])
        

        # WMAX = SC1 + SC2 (unit; m)
        self.thickUpp   = 0.3
        self.thickLow   = 1.2
        self.storCapUpp = self.thickUpp * \
                         (self.satVolWC1 - self.resVolWC1)
        self.storCapLow = self.thickLow * \
                         (self.satVolWC2 - self.resVolWC2)
        self.rootZoneWaterStorageCap = self.storCapUpp + self.storCapLow                      # This is called as WMAX in the original pcrcalc script. 

        # orographyBeta
        self.orographyBeta = vos.netcdf2PCRobjCloneWithoutTime(ncFile  = self.model_setup['topo_nc_file'],\ 
                                                               varName = "orographyBeta",\ 
                                                               cloneMapFileName  = self.cloneMap,\
                                                               LatitudeLongitude = True,\
                                                               specificFillValue = None,\
                                                               absolutePath = None)
        
        # land cover types
        self.coverTypes = ["forest", "grassland", "irrPaddy", "irrNonPaddy"]

        
        # fractions of natural land covers (forest and grassland) - this is over the entire cell area and not considering irrigated land
        fraction_forest    = vos.readPCRmapClone(v = self.model_setup["fraction_forest"], \
                                                      cloneMapFileName = self.clone, \
                                                      tmpDir           = self.tmp_folder, \
                                                      absolutePath     = None, \
                                                      isLddMap         = False, \
                                                      cover            = None, \
                                                      isNomMap         = False)
        fraction_grassland = vos.readPCRmapClone(v = self.model_setup["fraction_grassland"], \
                                                      cloneMapFileName = self.clone, \
                                                      tmpDir           = self.tmp_folder, \
                                                      absolutePath     = None, \
                                                      isLddMap         = False, \
                                                      cover            = None, \
                                                      isNomMap         = False)
        # correcting
        total_fraction_of_natural_before_correction = fraction_forest + fraction_grassland
        self.fraction_forest    = vos.getValDivZero(fraction_forest, total_fraction_of_natural_before_correction)
        self.fraction_grassland = pcr.max(0.0, 1.0 - self.fraction_forest)
        
        
        # irrTypeFracOverIrr = fraction each land cover type (paddy or nonPaddy) over the irrigation area (dimensionless) ; this value is constant for the entire simulation of the Aqueduct run
        fraction_paddy_over_irrigated_land     = vos.readPCRmapClone(v                = self.model_setup["fraction_paddy_over_irrigated_land"], \
                                                                     cloneMapFileName = self.clone, \
                                                                     tmpDir           = self.tmp_folder, \
                                                                     absolutePath     = None, \
                                                                     isLddMap         = False, \
                                                                     cover            = None, \
                                                                     isNomMap         = False)
        
        fraction_non_paddy_over_irrigated_land = vos.readPCRmapClone(v                = self.model_setup["fraction_non_paddy_over_irrigated_land"], \
                                                                     cloneMapFileName = self.clone, \
                                                                     tmpDir           = self.tmp_folder, \
                                                                     absolutePath     = None, \
                                                                     isLddMap         = False, \
                                                                     cover            = None, \
                                                                     isNomMap         = False)
        # - correcting
        total_irrigated_land_fraction_before_correction = fraction_paddy_over_irrigated_land + fraction_non_paddy_over_irrigated_land
        self.irrTypeFracOverIrr["irrPaddy"]    = vos.getValDivZero(fraction_paddy_over_irrigated_land, total_irrigated_land_fraction_before_correction)                                                            
        self.irrTypeFracOverIrr["irrNonPaddy"] = vos.getValDivZero(fraction_non_paddy_over_irrigated_land, total_irrigated_land_fraction_before_correction)                                                            


        

    def dynamic(self):

        # read the area/extent of irrigated lands
        self.dynamicIrrigationAreaFile = self.model_setup[" "]
        if self.dynamicIrrigationAreaFile.endswith(('.nc4','.nc')):
            fulldateInString = yearInString+"-01"+"-01"   
            self.irrigationArea = 10000. * pcr.cover(\
                 vos.netcdf2PCRobjClone(self.dynamicIrrigationAreaFile,\
                                            'irrigationArea',\
                     fulldateInString, useDoy = 'yearly',\
                             cloneMapFileName = self.cloneMap), 0.0)        # unit: m2 (input file is in hectare)
        
        # area of irrigation is limited by cellArea
        self.irrigationArea = pcr.max(self.irrigationArea, 0.0)              
        self.irrigationArea = pcr.min(self.irrigationArea, self.cell_area)  # limited by cellArea

        # calculate fracVegCover (for irrigation only): for "irrPaddy" and "irrNonPaddy"
        for coverType in self.coverTypes:
            if coverType.startswith('irr'):

                self.fractionArea[coverType] = 0.0    # reset 
                self.fractionArea[coverType] = self.irrTypeFracOverIrr[coverType]* self.irrigationArea # unit: m2
                self.fracVegCover[coverType] = pcr.min(1.0, self.fractionArea[coverType]/ self.cellArea) 

                # avoid small values
                self.fracVegCover[coverType] = pcr.rounddown(self.fracVegCover[coverType] * 1000.)/1000.


        # TODO: Please check the following!!
        
        # rescale land cover fractions (for all land cover types) - this is adopted from the function "scaleModifiedLandCoverFractions()" in the landSurface.py of PCR-GLOBWB
        # - calculate irrigatedAreaFrac (fraction of irrigation areas) 
        irrigatedAreaFrac = pcr.spatial(pcr.scalar(0.0))
        for coverType in self.coverTypes:
            if coverType.startswith('irr'):
                irrigatedAreaFrac = irrigatedAreaFrac + self.fracVegCover[coverType]

        totalArea  = pcr.spatial(pcr.scalar(0.0))
        totalArea += irrigatedAreaFrac

        # correction factor for forest and grassland (pristine Areas)
        lcFrac = pcr.max(0.0, 1.0 - totalArea)
        pristineAreaFrac = pcr.spatial(pcr.scalar(0.0))

        for coverType in self.coverTypes:         
            if not coverType.startswith('irr'):
                self.fracVegCover[coverType] = 0.0
                self.fracVegCover[coverType] = self.naturalFracVegCover[coverType] * lcFrac
                pristineAreaFrac             = pcr.cover(self.fracVegCover[coverType], 0.0)

        # check and make sure that totalArea = 1.0 for all cells
        totalArea += pristineAreaFrac
        totalArea = pcr.ifthen(self.landmask,totalArea)
        totalArea = pcr.cover(totalArea, 1.0)
        totalArea = pcr.ifthen(self.landmask,totalArea)
        a,b,c = vos.getMinMaxMean(totalArea - pcr.scalar(1.0))
        threshold = 1e-4
        if abs(a) > threshold or abs(b) > threshold:
            logger.error("fraction total (from all land cover types) is not equal to 1.0 ... Min %f Max %f Mean %f" %(a,b,c)) 


        irrigationArea
        
        # read land cover parameters
        land_cover_fraction = {}


        # read land cover parameters
        land_cover_fraction = {}
        minSoilDepthFrac = {}
        maxSoilDepthFrac = {}
        for coverType in coverTypes:
            
            # read land cover fractions
            land_cover_fraction[coverType] = 
            
            # read land cover parameters
            minSoilDepthFrac[coverType] =  
            maxSoilDepthFrac[coverType] = 
            
        # calculate tha aggregate land cover parameters
        minSoilDepthFrac_avg = pcr.scalar(0.0)
        maxSoilDepthFrac_avg = pcr.scalar(0.0)
        for coverType in coverTypes:
            minSoilDepthFrac_avg = minSoilDepthFrac_avg + minSoilDepthFrac[coverType] * land_cover_fraction[coverType]
            maxSoilDepthFrac_avg = maxSoilDepthFrac_avg + maxSoilDepthFrac[coverType] * land_cover_fraction[coverType]
            
        # arnoBeta
        self.arnoBeta = pcr.max(0.001,\
                 (maxSoilDepthFrac_avg-1.)/(1.-minSoilDepthFrac_avg)+\
                                           orographyBeta-0.01)                # Rens's line: BCF[TYPE]= max(0.001,(MAXFRAC[TYPE]-1)/(1-MINFRAC[TYPE])+B_ORO-0.01)

        # WMIN (unit: m): minimum local soil water capacity within the grid-cell
        self.rootZoneWaterStorageMin = minSoilDepthFrac_avg * self.rootZoneWaterStorageCap


        # re-calculate current model time using current pcraster timestep value
        self.modelTime.update(self.currentTimeStep())

            
        # calculating saturated area fraction
        if self.modelTime.isLastDayOfMonth():
            
            logger.info(" \n\n Calculating for time %s \n\n", self.modelTime.currTime)

            # read monthly S1 and S2 (m)
            monthly_storUpp    = vos.netcdf2PCRobjClone(self.model_setup["monthly_s1_file"], \
                                                                   "total_runoff", \
                                                                   str(self.modelTime.fulldate), \
                                                                   None, \
                                                                   self.clone)
            monthly_storLow    = vos.netcdf2PCRobjClone(self.model_setup["monthly_s2_file"], \
                                                                   "total_runoff", \
                                                                   str(self.modelTime.fulldate), \
                                                                   None, \
                                                                   self.clone)
            monthly_total_soil_storage = monthly_storUpp + monthly_storLow
            
            
            # calculate saturated area fraction
            saturated_area_fraction = 1.00 - \
             ((self.rootZoneWaterStorageCap - monthly_total_soil_storage) / (self.rootZoneWaterStorageCap - self.rootZoneWaterStorageMin))**(self.arnoBeta/(self.arnoBeta+1))
            
            
            # TODO: UNTIL THIS PART. We have to put all "self" things in the "initial" part.

            # reporting 
            # - time stamp for reporting
            timeStamp = datetime.datetime(self.modelTime.year,\
                                          self.modelTime.month,\
                                          self.modelTime.day,\
                                          0)
            logger.info("Reporting for time %s", self.modelTime.currTime)
            self.netcdf_report.data2NetCDF(self.model_setup["discharge_output_file"], \
                                           "discharge", \
                                           pcr.pcr2numpy(self.discharge, vos.MV), \
                                           timeStamp)


def main():

    model_setup = {}

    model_setup["clone_file"]     = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/lddsound_05min_version_20210330.map"
    
    model_setup["ldd_file"]       = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/lddsound_05min_version_20210330.map"
    
    model_setup["cell_area_file"] = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/cdo_gridarea_clone_global_05min_correct_lats.nc"

    model_setup['resVolWC1']      = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/vmcRes_average_1_global_05arcmin.nc"
    model_setup['resVolWC2']      = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/vmcRes_average_2_global_05arcmin.nc"
    model_setup['satVolWC1']      = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/vmcSat_average_1_global_05arcmin.nc"
    model_setup['satVolWC2']      = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/vmcSat_average_2_global_05arcmin.nc"

    model_setup['topo_nc_file']   = "/projects/0/dfguu/users/edwin/data/pcrglobwb_input_aqueduct/version_2021-09-16/general/topography_parameters_5min_april_2021_global_covered_with_zero.nc"
    
    model_setup["irrigationArea"] = "irrigated_areas_historical_1960-2019.nc" 
    
    # fractions of natural land covers (forest and grassland)
    model_setup["fraction_forest"]    = "general/vegf_tall.map" 
    model_setup["fraction_grassland"] = "general/vegf_short.map" 


    # fractions of paddy and non paddy over irrigated areas
    model_setup["fraction_paddy_over_irrigated_land"]     = "general/fractionPaddy_extrapolated.map" 
    model_setup["fraction_non_paddy_over_irrigated_land"] = "general/fractionNonPaddy_extrapolated.map" 
    





 dem_minimum dem_maximum dem_average dem_standard_deviation slopeLength orographyBeta tanslope dzRel0000 dzRel0001 dzRel0005 dzRel0010 dzRel0020 dzRel0030 dzRel0040 dzRel0050 dzRel0060 dzRel0070 dzRel0080 dzRel0090 dzRel0100
cdo    showname: Processed 20 variables [0.16s 81MB].



'resVolWC1',\
'resVolWC2',\
'satVolWC1',\
'satVolWC2'



    model_setup["monthly_s1_file"] = "/projects/0/managed_datasets/hypflowsci6_v1.0/output/gswp3-w5e5/historical-reference/pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_storUppTotal_global_monthly-average_1960_2019_basetier1.nc"
    model_setup["monthly_s2_file"] = "/projects/0/managed_datasets/hypflowsci6_v1.0/output/gswp3-w5e5/historical-reference/pcrglobwb_cmip6-isimip3-gswp3-w5e5_image-aqueduct_historical-reference_storLowTotal_global_monthly-average_1960_2019_basetier1.nc"
    
    model_setup["start_date"] = "1960-01-31"
    model_setup["end_date"]   = "2019-12-31"

    model_setup["output_dir"] = "/scratch/depfg/sutan101/discharge_30sec_gmd_paper/monthly_1958-2015_splitted/" + model_setup["start_date"] + "_to_" + model_setup["end_date"] + "/"

    model_setup["discharge_output_file"] = model_setup["output_dir"] + "/" + "discharge_30sec_monthAvg_" + model_setup["start_date"] + "_to_" + model_setup["end_date"] + ".nc"


    print(model_setup["output_dir"])
    
    # make output and temporary folders
    if os.path.exists(model_setup["output_dir"]): shutil.rmtree(model_setup["output_dir"])
    os.makedirs(model_setup["output_dir"])
    # - make temporary folder
    model_setup["tmp_dir"] = model_setup["output_dir"] +  "/tmp/"
    os.makedirs(model_setup["tmp_dir"])


    # logger
    # - making a log directory
    log_file_directory = model_setup["output_dir"] + "/" + "log/"
    os.makedirs(log_file_directory)
    # - initialize logging
    vos.initialize_logging(log_file_directory)
    

    # timeStep info: year, month, day, doy, hour, etc
    currTimeStep = ModelTime() 
    currTimeStep.getStartEndTimeSteps(model_setup["start_date"], model_setup["end_date"])
    

    # Running the deterministic_runner
    logger.info('Starting the calculation.')
    deterministic_runner = DeterministicRunner(currTimeStep, model_setup)
    dynamic_framework = DynamicFramework(deterministic_runner,currTimeStep.nrOfTimeSteps)
    dynamic_framework.setQuiet(True)
    dynamic_framework.run()  
    
        
if __name__ == '__main__':
    sys.exit(main())
