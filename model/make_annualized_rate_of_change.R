## #############################################################################
## 
## DETERMINE AROC, 2010-2017
## 
## Author: Nat Henry
## Created: August 22, 2019
## Purpose: Determine annualized rate of change from 2010 through 2017 by grid
##   cell (and aggregate) in India
## 
## #############################################################################

## SET INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

core_repo <- '{REDACTED}'
ig <- 'u5m'
rd <- '2019_06_07_rtr'
shp_version <- '2019_05_06'
# Age groups to analyze
age_groups <- c('under5','neonatal')
wp_pop_measure <- 'a0004t'
# GBD age group IDs associated with each age group
ag_ids <- list(
  under5=2:5, neonatal=2:3
)
# Years for which we have data
data_years <- 2000:2017
# Bounds for annualized rates of change
aroc_bounds <- list(
  c(2000, 2010),
  c(2000, 2017),
  c(2010, 2017)
)


## SETUP FOR ALL REGIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load MBG functions and packages
commondir <- paste0(core_repo, '/mbg_central/share_scripts/common_inputs/')
source(sprintf('%s/mbg_central/setup.R',core_repo))
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)
library(ncdf4, lib.loc='/home/j/temp/nathenry/u5m/visualization/inputs/packages/')
source('/ihme/cc_resources/libraries/current/r/get_age_metadata.R')
source('/ihme/cc_resources/libraries/current/r/get_location_metadata.R')

## Load admin metadata
ad2_meta <- as.data.table(foreign::read.dbf(
  get_admin_shapefile(admin_level=2, version=shp_version, suffix='.dbf')
))[, .(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME)]
ad1_meta <- unique(ad2_meta[, .(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME)])
ad0_meta <- unique(ad2_meta[, .(ADM0_CODE, ADM0_NAME)])


## Cell-level AROC calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(aroc_bound_group in aroc_bounds){
  aroc_y1 <- aroc_bound_group[1]
  aroc_y2 <- aroc_bound_group[2]
  message(glue::glue("\n\n** MAKING AROC FOR {aroc_y1} TO {aroc_y2} **\n"))

  ## Iterate by region
  for(region in regions){
    # Load simple raster
    gaul_list <- get_adm0_codes(region, shapefile_version=shp_version)
    simple_polygon_list <- load_simple_polygon(
      gaul_list = gaul_list, buffer = 1, tolerance = 0.4)
    subset_shape        <- simple_polygon_list[[1]]
    simple_polygon      <- simple_polygon_list[[2]]

    ## Load list of raster inputs (pop and simple)
    raster_list <- build_simple_raster_pop(subset_shape, link_table=shp_version)
    simple_raster <- raster_list[['simple_raster']]

    ## Iterate by age group
    for(ag in age_groups){
      indicator <- paste0('died_',ag)
      # Setup age group
      message(sprintf("Working on %s for %s...",region,indicator))
      sharedir <- sprintf('/share/geospatial/mbg/%s/%s/output/%s/',ig,indicator,rd)
      proj_dir <- paste0(sharedir,'/proj_ind/')
      dir.create(proj_dir, showWarnings=FALSE)

      # Load cell pred
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

      ## Keep first and last years of data
      if(aroc_y1 >= aroc_y2) stop("Second year must be after the first.")
      y1_idx <- which(aroc_y1==data_years)
      y2_idx <- which(aroc_y2==data_years)
      cp_start <- cp_reg[((ncell*(y1_idx-1))+1):(ncell*y1_idx),]
      cp_end <- cp_reg[((ncell*(y2_idx-1))+1):(ncell*y2_idx),]
      rm(cp_reg)
      gc(full=TRUE)
      ## Find annualized rate of change across the whole year range
      aroc_cp <- exp( log(cp_end/cp_start) / (aroc_y2-aroc_y1) )

      ## Create and save summary rasters
      message('  Summarizing AROC cell preds:')
      for(summeasure in c('mean','lower','upper')){
        message(sprintf('    %s',summeasure))
        tmp <-  suppressMessages(make_cell_pred_summary(
          draw_level_cell_pred = aroc_cp,
          mask = simple_raster,
          return_as_raster = TRUE,
          summary_stat = summeasure
        ))
        writeRaster(
          tmp-1, format='GTiff', overwrite=TRUE,
          file=sprintf(
            '%s/%s_%s_raked_%s_%s_%s_aroc.tif',proj_dir,indicator,summeasure,
            region, aroc_y1, aroc_y2
          )
        )
      } # End raster summarization loop
    } # End age group loop
  } # End region loop


  ## Admin-level AROC calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  aroc_calc <- function(val_start, val_end, tdiff){
    return(
      exp(log(val_end/val_start) / tdiff)
    )
  }


  message("Working on admin-level AROC calculations")
  ## Iterate by age group
  for(ag in age_groups){
    # Setup age group
    indicator <- paste0('died_',ag)
    sharedir <- sprintf('/share/geospatial/mbg/%s/%s/output/%s/',ig,indicator,rd)
    proj_dir <- paste0(sharedir,'/proj_ind/')
    dir.create(proj_dir, showWarnings=FALSE)

    ad_list <- list(ad0=list(), ad1=list(), ad2=list())
    for(region in regions){
      message(sprintf("Working on %s for %s...",region,indicator))

      ## Iterate by admin level
      for(ad_lev in 0:2){
        message(sprintf('  Summarizing Admin%i AROC...', ad_lev))
        # Load annual draws
        yr_draws <- fread(sprintf(
          '%s/%s_%s_draws_raked_ad%s.csv',sharedir,indicator,region,ad_lev
        ))
        if(names(yr_draws)[1]=='V1') yr_draws[,1 := NULL]
        draw_cols <- grep('^V',names(yr_draws),value=TRUE)
        ndraws <- length(draw_cols)
        # Keep years involved in AROC calculation; merge and get AROC
        keep_cols <- c('name',draw_cols)
        draws_y1 <- yr_draws[year==aroc_y1, ..keep_cols]
        draws_y2 <- yr_draws[year==aroc_y2, ..keep_cols]
        draws_merged <- merge(
          x=draws_y1, y=draws_y2, by='name', all=T, suffixes=c('_y1','_y2')
        )
        aroc_cols <- paste0('aroc_',1:ndraws)
        tdiff <- aroc_y2 - aroc_y1
        for(idx in 1:ndraws){
          draws_merged[[aroc_cols[idx]]] <- aroc_calc(
            val_start = draws_merged[[sprintf('V%i_y1',idx)]],
            val_end = draws_merged[[sprintf('V%i_y2',idx)]],
            tdiff = tdiff
          )
        }
        keep_cols <- c('name',aroc_cols)
        draws_merged <- draws_merged[, ..keep_cols]
        draw_summary <- data.table(name=draws_merged$name)
        draw_summary$mean <- rowMeans(draws_merged[, ..aroc_cols], na.rm=T)-1
        draw_summary$lower <- apply(
          draws_merged[, ..aroc_cols], 1, function(x) quantile(x,0.025,na.rm=T)-1
        )
        draw_summary$upper <- apply(
          draws_merged[, ..aroc_cols], 1, function(x) quantile(x,0.975,na.rm=T)-1
        )
        setnames(draw_summary, 'name', sprintf('ADM%i_CODE',ad_lev))

        ## Add summary data.table to full list!
        ad_list[[paste0('ad',ad_lev)]][[region]] <- draw_summary
      } # END admin level loop
    } # End region loop

    ## Combine summary lists into full data.tables and save!
    ad0_summ_full <- rbindlist(ad_list$ad0)
    ad1_summ_full <- rbindlist(ad_list$ad1)
    ad2_summ_full <- rbindlist(ad_list$ad2)
    proj_out_templ <- sprintf(
      '%s/%s_raked_aroc_ad%%i_fullsummary.csv',proj_dir,indicator
    )
    write.csv(ad0_summ_full, file=sprintf(proj_out_templ,0), row.names=FALSE)
    write.csv(ad1_summ_full, file=sprintf(proj_out_templ,1), row.names=FALSE)
    write.csv(ad2_summ_full, file=sprintf(proj_out_templ,2), row.names=FALSE)

    # Add on metadata; save age group subset
    ad0_ind <- merge(x=ad0_meta[ADM0_NAME=='India',],y=ad0_summ_full,by='ADM0_CODE')
    ad1_ind <- merge(x=ad1_meta[ADM0_NAME=='India',],y=ad1_summ_full,by='ADM1_CODE')
    ad2_ind <- merge(x=ad2_meta[ADM0_NAME=='India',],y=ad2_summ_full,by='ADM2_CODE')
    collab_out_dir <- '{REDACTED}'
    clb_tmpl <- '{REDACTED}'
    write.csv(ad0_ind, file=sprintf(clb_tmpl,0), row.names=FALSE)
    write.csv(ad1_ind, file=sprintf(clb_tmpl,1), row.names=FALSE)
    write.csv(ad2_ind, file=sprintf(clb_tmpl,2), row.names=FALSE)
  } # End age group loop


  ## SUMMARIZE RASTER PROJECTIONS BY REGION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for(ag in age_groups){
    indicator <- paste0('died_',ag)
    proj_dir <- '{REDACTED}'
    ## CREATE RASTER SUMMARIES
    for(summ_measure in c('mean','lower','upper')){
      rast_fps <- list.files(
        proj_dir, ignore.case=TRUE, full.names=TRUE,
        pattern=sprintf(
          '%s_%s_raked_[a-z]+[\\+]?[a-z]*_%s_%s_aroc.tif',
          indicator,summ_measure,aroc_y1,aroc_y2
        )
      )
      rast_combined <- do.call('merge', lapply(rast_fps, function(x) raster::brick(x)))
      writeRaster(
        rast_combined, format='GTiff', overwrite=TRUE,
        file=sprintf(
          '%s/%s_%s_raked_%s_%s_aroc.tif',
          proj_dir,indicator,summ_measure,aroc_y1,aroc_y2
        )
      )
    } # End summary measure loop
  } # End age group loop

  message(glue::glue("... Saved aroc for {aroc_y1}-{aroc_y2} successfully."))

} # End full loop for year set


message("~~~~ FIN ~~~~")