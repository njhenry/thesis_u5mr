## #######################################################################################
##
## MATCH SRS DATA TO DISTRICT GROUPINGS
## AUTHOR: Nat Henry
## PURPOSE: SRS data is reported by NSSO "natural division," which align with groups of
##   districts within a single state. This script matches the SRS 2017 results for infant
##   mortality, reported by "natural division," with their component districts
##
## #######################################################################################

## SETUP -------------------------------------------------------------------------------->

library(data.table); library(foreign); library(knitr)

work_dir <- '{REDACTED}'
shp_dir <- '{REDACTED}'

# Load natural divisions
nds <- data.table::fread(file.path(work_dir, 'natural_divisions_shp_matched.csv'))
# Load SRS IMR data
imr_raw <- data.table::fread(file.path(work_dir, 'srs_imr_by_nd_2017.csv'))
# Load admin2 shapefile metadata
shp_meta <- data.table::as.data.table(foreign::read.dbf(
  file.path(shp_dir,'IND_full_ad2_for_mapping.dbf')
))

## Part 1: Match IMR data to NSSO table ------------------------------------------------->

imr_raw[, nd_name_imr := nd_name ]
nds[, nd_name_index := nd_name ]
test <- merge(
  x=imr_raw[, .(ADM1_NAME, nd_code, nd_name_imr)],
  y=unique(nds[, .(ADM1_NAME, nd_code, nd_name_index)]),
  all=TRUE
)
test[is.na(nd_name_index), ]
test[!is.na(nd_name_imr) & (nd_name_index != nd_name_imr), ]

keep_states <- test[!is.na(nd_name_imr), unique(ADM1_NAME)]
test[ADM1_NAME %in% keep_states & is.na(nd_name_imr), ]

## Part 2: Merge district IDs from admin2 shapefile onto NSSO table --------------------->

nsso_shp_match <- merge(
  x = shp_meta[, .(ADM1_NAME, ADM2_NAME, DIST_ID)][, shp := 1],
  y = nds[, .(ADM1_NAME,district_name)][, nsso := 1],
  by.x=c("ADM1_NAME", "ADM2_NAME"), by.y=c("ADM1_NAME", "district_name"), all=T
)[, ADM1_NAME := as.character(ADM1_NAME) ][order(ADM1_NAME, ADM2_NAME)]
nsso_check <- nsso_shp_match[ (is.na(nsso) | is.na(shp)), ]

# Check only states where IMR was reported by region
imr_states <- unique(imr_raw$ADM1_NAME)
imr_states[imr_states == 'Jammu & Kashmir'] <- 'Jammu and Kashmir'
imr_states[imr_states == 'Orissa'] <- 'Odisha'
imr_states[imr_states == 'Chattisgarh'] <- 'Chhattisgarh'

for(imr_state in imr_states){
  nsso_check_sub <- nsso_check[ ADM1_NAME == imr_state, ]
  if(nrow(nsso_check_sub)){
    print(knitr::kable(nsso_check_sub))
    message('\n\n\n')
  }
}

setnames(shp_meta, 'ADM2_NAME', 'district_name')
nds[shp_meta, DIST_ID := i.DIST_ID, on=c('ADM1_NAME','district_name')]
data.table::fwrite(nds, file=file.path(work_dir, 'natural_divisions_shp_matched.csv'))
