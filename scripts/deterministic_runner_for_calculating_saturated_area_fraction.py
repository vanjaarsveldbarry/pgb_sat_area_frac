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

import ncConverter as netcdf_writer
import virtualOS as vos

import logging
import configparser
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

        # set the landmask based on the extent of ldd
        self.landmask = pcr.ifthen(pcr.defined(self.ldd), pcr.boolean(1.0))

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
        self.netcdf_report.createNetCDF(self.model_setup["saturated_area_fraction_output_file"],\
                                        "satAreaFrac",\
                                        "-")
        
    def initial(self): 
        
        # read soil parameters from the netCDF file:
        soilParameters = ['resVolWC1',\
                          'resVolWC2',\
                          'satVolWC1',\
                          'satVolWC2']
        for var in soilParameters:
            soil_input_file = self.model_setup[var]
            vars(self)[var] = vos.readPCRmapClone(soil_input_file, self.clone, self.tmp_folder)
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
                                                               cloneMapFileName  = self.clone,\
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
        self.fraction_forest     = vos.getValDivZero(fraction_forest, total_fraction_of_natural_before_correction)
        self.fraction_grassland  = pcr.max(0.0, 1.0 - self.fraction_forest)
        self.naturalFracVegCover = {}
        self.naturalFracVegCover["forest"]    = self.fraction_forest
        self.naturalFracVegCover["grassland"] = self.fraction_grassland
        
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
        self.irrTypeFracOverIrr = {}
        self.irrTypeFracOverIrr["irrPaddy"]    = vos.getValDivZero(fraction_paddy_over_irrigated_land, total_irrigated_land_fraction_before_correction)                                                            
        self.irrTypeFracOverIrr["irrNonPaddy"] = vos.getValDivZero(fraction_non_paddy_over_irrigated_land, total_irrigated_land_fraction_before_correction)                                                            
        
        # read some land cover parameters, related to Arno scheme
        self.minSoilDepthFrac = {}  
        self.maxSoilDepthFrac = {}
        for coverType in self.coverTypes:
            self.minSoilDepthFrac[coverType] = vos.readPCRmapClone(v = self.model_setup["minSoilDepthFrac_" + coverType], \
                                                                   cloneMapFileName = self.clone, \
                                                                   tmpDir           = self.tmp_folder, \
                                                                   absolutePath     = None, \
                                                                   isLddMap         = False, \
                                                                   cover            = None, \
                                                                   isNomMap         = False)  
            self.maxSoilDepthFrac[coverType] = vos.readPCRmapClone(v = self.model_setup["maxSoilDepthFrac_" + coverType], \
                                                                   cloneMapFileName = self.clone, \
                                                                   tmpDir           = self.tmp_folder, \
                                                                   absolutePath     = None, \
                                                                   isLddMap         = False, \
                                                                   cover            = None, \
                                                                   isNomMap         = False)

    def dynamic(self):

        coverTypes = ["forest", "grassland", "irrPaddy", "irrNonPaddy"]
        
        # re-calculate current model time using current pcraster timestep value
        self.modelTime.update(self.currentTimeStep())

        # at the beginning of every year, read/set land cover parameters (related to the Arno Scheme) based on the exents of irrigation areas and natural land cover areas
        if self.modelTime.isFirstDayOfMonth() or self.modelTime.isFirstTimestep():

            # read the area/extent of irrigated lands
            self.dynamicIrrigationAreaFile = self.model_setup["irrigationArea"]
            if self.dynamicIrrigationAreaFile.endswith(('.nc4','.nc')):
                fulldateInString = str(self.modelTime.year) + "-01"+"-01"   
                self.irrigationArea = 10000. * pcr.cover(\
                     vos.netcdf2PCRobjClone(self.dynamicIrrigationAreaFile,\
                                                'irrigationArea',\
                         fulldateInString, useDoy = 'yearly',\
                                 cloneMapFileName = self.clone), 0.0)        # unit: m2 (input file is in hectare)
            
            # area of irrigation is limited by cellArea
            self.irrigationArea = pcr.max(self.irrigationArea, 0.0)              
            self.irrigationArea = pcr.min(self.irrigationArea, self.cell_area)  # limited by cellArea
		    
            # calculate fracVegCover (for irrigation only): for "irrPaddy" and "irrNonPaddy"
            self.fractionArea = {}
            self.fracVegCover = {}
            for coverType in self.coverTypes:
                if coverType.startswith('irr'):
                    self.fractionArea[coverType] = self.irrTypeFracOverIrr[coverType]* self.irrigationArea # unit: m2
                    self.fracVegCover[coverType] = pcr.min(1.0, self.fractionArea[coverType] / self.cell_area) 
                    # avoid small values
                    self.fracVegCover[coverType] = pcr.rounddown(self.fracVegCover[coverType] * 1000.)/1000.
            
            # rescale land cover fractions (for all land cover types)
            irrigatedAreaFrac = pcr.spatial(pcr.scalar(0.0))
            for coverType in self.coverTypes:
                if coverType.startswith('irr'):
                    irrigatedAreaFrac = irrigatedAreaFrac + self.fracVegCover[coverType]
		    
            # total area fraction after irrigatedAreaFrac
            totalArea  = pcr.spatial(pcr.scalar(0.0))
            totalArea += irrigatedAreaFrac
		    
            # natural/pristine area fraction (for forest and grassland)
            lcFrac = pcr.max(0.0, 1.0 - totalArea)
            pristineAreaFrac = pcr.spatial(pcr.scalar(0.0))
            for coverType in self.coverTypes:         
                if not coverType.startswith('irr'):
                    self.fracVegCover[coverType] = self.naturalFracVegCover[coverType] * lcFrac
                    pristineAreaFrac += self.fracVegCover[coverType]
		    
            # Add pristine area fractions to the total
            totalArea += pristineAreaFrac

            # Ensure totalArea is exactly 1.0 across all cells
            totalArea = pcr.max(0.0, pcr.min(1.0, totalArea))

            # Optionally, rescale land cover fractions so that they sum exactly to 1.0
            scale_factor = pcr.ifthen(self.landmask, pcr.scalar(1.0) / totalArea)
            for coverType in self.coverTypes:
                self.fracVegCover[coverType] = self.fracVegCover[coverType] * scale_factor

            # Check if totalArea equals 1
            totalArea = pcr.ifthen(self.landmask, totalArea)
            a, b, c = vos.getMinMaxMean(totalArea - pcr.scalar(1.0))
            threshold = 1e-4
            if abs(a) > threshold or abs(b) > threshold:
                logger.error("fraction total (from all land cover types) is not equal to 1.0 ... Min %f Max %f Mean %f" % (a,b,c)) 
            
        # Calculating saturated area fraction (the formula can be updated here)
        if self.modelTime.isLastDayOfMonth():
            
            logger.info(" \n\n Calculating for time %s \n\n", self.modelTime.currTime)

            # read monthly S1 and S2 (m)
            monthly_storUpp = vos.netcdf2PCRobjClone(self.model_setup["monthly_s1_file"], "upper_soil_storage", str(self.modelTime.fulldate), None, self.clone)
            monthly_storLow = vos.netcdf2PCRobjClone(self.model_setup["monthly_s2_file"], "lower_soil_storage", str(self.modelTime.fulldate), None, self.clone)
            monthly_total_soil_storage = monthly_storUpp + monthly_storLow
            
            # Calculate saturated area fraction with the provided formula
            saturated_area_fraction = (monthly_storUpp + monthly_storLow) / (self.storCapUpp + self.storCapLow)
            
            # Reporting 
            timeStamp = datetime.datetime(self.modelTime.year, self.modelTime.month, self.modelTime.day, 0)
            logger.info("Reporting for time %s", self.modelTime.currTime)
            self.netcdf_report.data2NetCDF(self.model_setup["saturated_area_fraction_output_file"], "satAreaFrac", pcr.pcr2numpy(saturated_area_fraction, vos.MV), timeStamp)


def main():

    model_setup = {}
    # Read the ini file
    config = configparser.ConfigParser()
    config.read(sys.argv[1])

    model_setup["output_dir"] = sys.argv[2]
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

    model_setup["input_dir"] = config.get('globalOptions', 'inputDir')
    model_setup["clone_file"] = f"{model_setup['input_dir']}/{config.get('globalOptions', 'cloneMap')}"
    model_setup["ldd_file"] = f"{model_setup['input_dir']}/{config.get('globalOptions', 'ldd_file')}"
    model_setup["cell_area_file"] = f"{model_setup['input_dir']}/{config.get('globalOptions', 'cell_area_file')}"
    model_setup['resVolWC1']      = f"{model_setup['input_dir']}general/vmcRes_average_1_global_05arcmin.nc"
    model_setup['resVolWC2']      = f"{model_setup['input_dir']}general/vmcRes_average_2_global_05arcmin.nc"
    model_setup['satVolWC1']      = f"{model_setup['input_dir']}general/vmcSat_average_1_global_05arcmin.nc"
    model_setup['satVolWC2']      = f"{model_setup['input_dir']}general/vmcSat_average_2_global_05arcmin.nc"
    model_setup['topo_nc_file']   = f"{model_setup['input_dir']}general/topography_parameters_5min_april_2021_global_covered_with_zero.nc"
    model_setup["irrigationArea"] = f"{model_setup['input_dir']}historical_and_ssp_files/irrigated_areas_historical_1960-2019.nc" 
    model_setup["fraction_forest"]    = f"{model_setup['input_dir']}general/vegf_tall.map" 
    model_setup["fraction_grassland"] = f"{model_setup['input_dir']}general/vegf_short.map" 
    model_setup["fraction_paddy_over_irrigated_land"]     = f"{model_setup['input_dir']}general/fractionPaddy_extrapolated.map" 
    model_setup["fraction_non_paddy_over_irrigated_land"] = f"{model_setup['input_dir']}general/fractionNonPaddy_extrapolated.map" 
    model_setup["minSoilDepthFrac_forest"]      = f"{model_setup['input_dir']}general/minf_tall.map"
    model_setup["minSoilDepthFrac_grassland"]   = f"{model_setup['input_dir']}general/minf_short.map"
    model_setup["minSoilDepthFrac_irrPaddy"]    = f"0.99"
    model_setup["minSoilDepthFrac_irrNonPaddy"] = f"0.99"
    model_setup["maxSoilDepthFrac_forest"]      = f"{model_setup['input_dir']}general/maxf_tall.map"
    model_setup["maxSoilDepthFrac_grassland"]   = f"{model_setup['input_dir']}general/maxf_short.map"
    model_setup["maxSoilDepthFrac_irrPaddy"]    = f"1.01"
    model_setup["maxSoilDepthFrac_irrNonPaddy"] = f"1.01"

    model_setup["start_date"] = sys.argv[3]
    model_setup["end_date"]   = sys.argv[4]
    model_setup["monthly_s1_file"] = sys.argv[5]
    model_setup["monthly_s2_file"] = sys.argv[6]

    model_setup["saturated_area_fraction_output_file"] = model_setup["output_dir"] + f"sat_area_fraction.nc"

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
