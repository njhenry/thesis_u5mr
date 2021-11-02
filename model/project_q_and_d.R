## #############################################################################
## 
## PROJECT Q AND D IN INDIA
## 
## Author: Nat Henry
## Created: August 1, 2019
## Purpose: Project information for q (mortality probability) and d (deaths)
##   to 2025 and 2030 for the India mortality collaboration
## 
## #############################################################################

## SET INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

in_dir <- '{REDACTED}'
core_repo <- '{REDACTED}'
ig <- 'u5m'
rd <- '2019_06_07_rtr'
shp_version <- '2019_05_06'
# Age groups to analyze
age_groups <- c('under5','neonatal')
wp_pop_measure <- 'a0004t'
ag_ids <- list(
  under5=2:5, infant=2:4, neonatal=2:3
)
data_years <- 2000:2017
proj_years <- c(2025, 2030)

## SETUP FOR ALL REGIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load MBG functions and packages
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)
library(ncdf4, lib.loc='{REDACTED}')
source('{REDACTED}')
source('{REDACTED}')


## PROJECT CELL PRED + CREATE UNRAKED ADMIN ESTIMATES ~~~~~~~~~~~~~~~~~~~~~~~~~~

## Iterate by region
for(region in regions){
  # Load simple raster
  simple_raster_fp <- '{REDACTED}'
  attach(simple_raster_fp) # Create a temporary environment from this Rdata object
  simple_raster <- simple_raster # Add an object from the temporary env to your working env
  detach()

  ## Iterate by age group
    for(ag in age_groups){
    indicator <- paste0('died_',ag)
    # Setup age group
    message(sprintf("Working on %s for %s...",region,indicator))
    sharedir <- '{REDACTED}'
    proj_dir <- paste0(sharedir,'/proj_ind/')
    dir.create(proj_dir, showWarnings=FALSE)
    # Create dummy "raking factors" data.table with RFs of 1
    ad0_codes <- get_adm0_codes(region)
    gbd_loc_ids <- load_adm0_lookup_table()[gadm_geoid %in% ad0_codes, loc_id]
    rf_table <- CJ(
      location_id = gbd_loc_ids,
      year = proj_years,
      rf = 1
    )
    write.csv(
      rf_table, row.names=FALSE,
      file=sprintf('%s/%s_%s_rf.csv',proj_dir,indicator,region)
    )
    pop_scalar_table <- CJ(
      location_id = gbd_loc_ids,
      year = proj_years,
      pop_scalar = 1
    )
    write.csv(
      pop_scalar_table, row.names=FALSE,
      file=sprintf('%s/%s_%s_pop_rf.csv',proj_dir,indicator,region)
    )

    cp_reg <- readRDS(sprintf(
      '%s/%s_raked_cell_draws_eb_bin0_%s_0.RDS',sharedir,indicator,region
    ))
    ## Check cell pred size
    ndraw <- dim(cp_reg)[2]
    ncell <- length(cellIdx(simple_raster))
    if(dim(cp_reg)[1] != ncell * length(data_years)){
      if(dim(cp_reg)[1] == (ncell+1)*length(data_years)){
        ## FIX YOUR SIMPLE RASTER
        message("Simple raster off by one... finding extra NAs and fixing")
        ## Load cell pred and keep first and last years
        if(not('cp_reg_ur' %in% ls())){
          cp_reg_ur <- readRDS(sprintf(
            '%s/%s_cell_draws_eb_bin0_%s_0.RDS',sharedir,indicator,region
          ))
        }
        if(nrow(cp_reg)-nrow(cp_reg_ur)!=length(data_years)) stop("Fix didn't work.")
        find_missing <- data.table(
          unraked = is.na(cp_reg_ur[1:(ncell+1),1]),
          raked = is.na(cp_reg[1:(ncell+1),1])
        )
        find_missing[, id := .I ]
        added_row <- first(find_missing[unraked==0 & raked==1,id])
        ## Check that rows across all years have this same issue
        for(i in 1:(length(data_years)-1)){
          if(not(
            !is.na(cp_reg_ur[(ncell*i)+added_row,1]) & 
            is.na(cp_reg[((ncell+1)*i)+added_row,1])
          )) stop("ISSUE")
        }
        # Drop the missing index from the raked cell pred
        drop_rows <- 0:(length(data_years)-1)*(ncell+1)+added_row
        if(not(all(is.na(cp_reg[drop_rows,])))) stop("Don't delete non-NA data")
        cp_reg <- cp_reg[-drop_rows,]
        if(dim(cp_reg)[1] != ncell * length(data_years)) stop("Fix failed.")
        message("Fix worked! Cleaning up...")
        rm('cp_reg_ur')
        gc(full=TRUE)
      } else {
        stop("Issue with cell pred dimensions")
      }
    }
    ## Keep first and last years
    cp_ly <- cp_reg[(ncell*(length(data_years)-1)+1):nrow(cp_reg),]
    cp_fy <- cp_reg[1:ncell,]
    rm(cp_reg)
    gc(full=TRUE)
    ## Find annualized rate of change across the whole year range
    aroc_cp <- exp( log(cp_ly/cp_fy) / diff(range(data_years)) )
    ## Project to future years
    proj_cp <- do.call(
      rbind,
      lapply( proj_years, function(yr) cp_ly * (aroc_cp ^ (yr - max(data_years))) )
    )

    ## Create and save summary rasters
    message('  Summarizing projected cell preds:')
    for(summeasure in c('mean','lower','upper')){
      message(sprintf('    %s',summeasure))
      tmp <-  suppressMessages(make_cell_pred_summary(
        draw_level_cell_pred = proj_cp,
        mask = simple_raster,
        return_as_raster = TRUE,
        summary_stat = summeasure
      ))
      writeRaster(
        tmp, format='GTiff', overwrite=TRUE,
        file=sprintf(
          '%s/%s_%s_unraked_%s_%s_%s_projections.tif',proj_dir,indicator,summeasure,
          region, min(proj_years),max(proj_years)
        )
      )
    }

    ## Aggregate to admin2/admin1/admin0
    pop_release <- '2017_04_27'
    message("  Fractionally aggregating to admin2/1/0 and saving")
    fractional_agg_rates(
      cell_pred = proj_cp,
      simple_raster = simple_raster,
      simple_polygon = simple_raster,
      pixel_id = cellIdx(simple_raster),
      shapefile_version = shp_version,
      reg = region,
      pop_measure = wp_pop_measure,
      year_list = proj_years,
      use_intermediate_years = FALSE,
      interval_mo = 12,
      rake_subnational = FALSE,
      sharedir = sprintf('/share/geospatial/mbg/u5m/%s/',indicator),
      run_date = rd,
      indicator = indicator,
      main_dir = sharedir,
      rake_method = "linear",
      age = 0,
      holdout = 0,
      return_objects = FALSE,
      countries_not_to_subnat_rake = countries_not_to_subnat_rake,
      custom_output_folder = proj_dir
    )
    rm('tmp','proj_cp','aroc_cp','cp_fy','cp_ly')
    gc(full=TRUE)
  }
}


## SUMMARIZE PROJECTIONS BY REGION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(ag in age_groups){
  ## CREATE ADMIN SUMMARIES
  indicator <- paste0('died_',ag)
  proj_dir <- '{REDACTED}'
  ad0_summ_list <- vector('list', length=length(regions))
  ad1_summ_list <- vector('list', length=length(regions))
  ad2_summ_list <- vector('list', length=length(regions))

  for(region in regions){
    message(sprintf("Summarizing %s for %s...",region,indicator))
    # Load admin data
    load(sprintf(
      '%s/%s_raked_admin_draws_eb_bin0_%s_0.RData',proj_dir,indicator,region
    ))
    ## Iterate through admin levels
    # Admin0
    draw_cols <- grep('^V',names(admin_0),value=TRUE)
    ad0_draws <- as.matrix(admin_0[, ..draw_cols])
    ad0_summary <- admin_0[, .(ADM0_CODE,year)]
    ad0_summary$mean <- rowMeans(ad0_draws, na.rm=TRUE)
    ad0_summary$lower <- apply(ad0_draws,1,function(x) quantile(x,0.025,na.rm=TRUE))
    ad0_summary$upper <- apply(ad0_draws,1,function(x) quantile(x,0.975,na.rm=TRUE))
    # Admin1
    draw_cols <- grep('^V',names(admin_1),value=TRUE)
    ad1_draws <- as.matrix(admin_1[, ..draw_cols])
    ad1_summary <- admin_1[, .(ADM1_CODE,year)]
    ad1_summary$mean <- rowMeans(ad1_draws,na.rm=TRUE)
    ad1_summary$lower <- apply(ad1_draws,1,function(x) quantile(x,0.025,na.rm=TRUE))
    ad1_summary$upper <- apply(ad1_draws,1,function(x) quantile(x,0.975,na.rm=TRUE))
    # Admin2
    draw_cols <- grep('^V',names(admin_2),value=TRUE)
    ad2_draws <- as.matrix(admin_2[, ..draw_cols])
    ad2_summary <- admin_2[, .(ADM2_CODE,year)]
    ad2_summary$mean <- rowMeans(ad2_draws,na.rm=TRUE)
    ad2_summary$lower <- apply(ad2_draws,1,function(x) quantile(x,0.025,na.rm=TRUE))
    ad2_summary$upper <- apply(ad2_draws,1,function(x) quantile(x,0.975,na.rm=TRUE))
    # Add to summary tables
    ad0_summ_list[[region]] <- ad0_summary
    ad1_summ_list[[region]] <- ad1_summary
    ad2_summ_list[[region]] <- ad2_summary
  }
  write.csv(
    rbindlist(ad0_summ_list), row.names=FALSE,
    file=sprintf('%s/%s_unraked_%s_%s_ad0_fullsummary.csv',proj_dir,indicator,min(proj_years),max(proj_years))
  )
  write.csv(
    rbindlist(ad1_summ_list), row.names=FALSE,
    file=sprintf('%s/%s_unraked_%s_%s_ad1_fullsummary.csv',proj_dir,indicator,min(proj_years),max(proj_years))
  )
  write.csv(
    rbindlist(ad2_summ_list), row.names=FALSE,
    file=sprintf('%s/%s_unraked_%s_%s_ad2_fullsummary.csv',proj_dir,indicator,min(proj_years),max(proj_years))
  )

  ## CREATE RASTER SUMMARIES
  for(summ_measure in c('mean','lower','upper')){
    rast_fps <- list.files(
      proj_dir, ignore.case=TRUE, full.names=TRUE,
      pattern=sprintf(
        '%s_%s_unraked_[a-z]+[\\+]?[a-z]*_%s_%s_projections.tif',
        indicator,summ_measure,min(proj_years),max(proj_years)
      )
    )
    rast_combined <- do.call('merge', lapply(rast_fps, function(x) raster::brick(x)))
    writeRaster(
      rast_combined, format='GTiff', overwrite=TRUE,
      file=sprintf(
        '%s/%s_%s_unraked_%s_%s_projections.tif',
        proj_dir,indicator,summ_measure,min(proj_years),max(proj_years)
      )
    )
  }
}
