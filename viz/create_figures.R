## #######################################################################################
##
## Create figures for thesis chapter on under-5 mortality in India
## CREATED: 26 July 2021
## AUTHOR: Nat Henry, github:@njhenry
##
## #######################################################################################

## SETUP -------------------------------------------------------------------------------->

# Load libraries
required_libs <- c('data.table','dplyr','ggplot2','glue','grid','gridExtra','RColorBrewer','scales','sf')
invisible(lapply(required_libs, library, character.only=TRUE))

# Set file paths
in_dir <- '{REDACTED}'
out_dir <- file.path(in_dir, 'pubs')
dir.create(out_dir, showWarnings=FALSE)

# Load mortality data
age_groups <- c('neonatal','infant','under5')
read_ags <- function(templ) fread(file.path(in_dir, templ))
summ <- lapply(age_groups, function(ag) read_ags(glue("india_{ag}_summary_complete.csv")))
proj <- lapply(age_groups, function(ag) read_ags(glue("died_{ag}_india_ad2_raked_2025_2030_20190826.csv")))
aroc <- lapply(
  age_groups[age_groups != 'infant'],
  function(ag) read_ags(glue("died_{ag}_raked_aroc_2000_2017_ad2_fullsummary.csv"))
)
names(summ) <- names(proj) <- age_groups
names(aroc) <- age_groups[age_groups != 'infant']

# Load admin0, admin1, and admin2 shapefiles
shp_dir <- '{REDACTED}'
shps <- list(
  ad1 = sf::st_read(file.path(shp_dir, "IND_full_ad1_for_mapping_no_islands.shp")),
  ad2 = sf::st_read(file.path(shp_dir, "IND_full_ad2_for_mapping.shp"))
)
admin_levels <- names(shps)

# Load SRS data by natural division and index linking natural divisions to districts
srs_tabs <- fread(file.path(in_dir, 'srs_2017', 'srs_imr_by_nd_2017.csv'))
nd_dist <- fread(file.path(in_dir, 'srs_2017', 'natural_divisions_shp_matched.csv'))
# Create a unique identifier for natural divisions
nd_idx <- unique(nd_dist[, .(ADM1_CODE, nd_code)])[order(ADM1_CODE, nd_code)][, nd_idx := .I]
nd_dist[nd_idx, nd_idx := i.nd_idx, on=c('ADM1_CODE', 'nd_code')]

# Construct the 'natural division' shapefile
nd_sf <- merge(x=shps$ad2[, c('DIST_ID')], y=nd_dist, on='DIST_ID')
nd_sf <- aggregate(
  nd_sf[, c('nd_idx')], by=list(nd_sf$nd_idx),
  FUN=first, simplify=FALSE
)
nd_sf$nd_idx <- unlist(nd_sf$nd_idx)
nd_sf <- merge(
  x=nd_sf,
  y=unique(nd_dist[, .(ADM1_CODE, ADM1_NAME, nd_name, nd_code, nd_idx)])
)
nd_sf_sub <- merge(x=nd_sf, y=srs_tabs[, c('ADM1_CODE', 'nd_code')])



## MAKE PLOTS --------------------------------------------------------------------------->

## ** (FIG 1: Maps of U5MR, IMR, and NMR across India in 2000 and 2017) **

color_scheme <- c(
  "#2D4386", "#0484ae", "#72ca88", "#ffecb3", "#F9BC87", "#f28c5a", "#d16161", "#862d2d"
)
breaks <- list(under5 = c(0,23,50,75), infant = c(0,20,40,60), neonatal = c(0,16,30,50))
titles <- list(
  neonatal = 'Neonatal (<1 month)', infant = 'Infant (<1 year)', under5 = 'Child (<5 years)'
)
india_refsys <- sf::st_crs(7755)
map_linecol = '#222222'
map_res <- 200

# Function to make an admin map by age group
make_admin_map <- function(ag){
  map_data <- summ[[ag]][(admin_level==2) & (year %in% c(2000, 2017)), ]
  map_labs <- breaks[[ag]]
  map_breaks <- map_labs / 1000
  col_max <- max(map_breaks)
  col_min <- min(map_breaks)
  map_data[mean_q > col_max, mean_q := col_max - 1E-5 ]
  map_data[mean_q < col_min, mean_q := col_min + 1E-5 ]
  plot_shp <- shps$ad2
  plot_shp$ADM2_CODE <- NULL
  data_merged <- merge(
    x = shps$ad2, y = map_data[, .(ADM2_CODE, year, mean_q)],
    by.x = 'DIST_ID', by.y = 'ADM2_CODE'
  )
  # Plot it
  legend_ag <- list(neonatal='Neonatal',infant='Infant',under5='Under-5')[[ag]]
  fig <- ggplot() +
    facet_wrap('year', ncol=2) +
    geom_sf(data = data_merged, aes(fill=mean_q), lwd=0.03, color=map_linecol) +
    geom_sf(data = shps$ad1, lwd=0.2, fill=NA, color=map_linecol) +
    scale_fill_gradientn(
      colors = color_scheme, breaks = map_breaks, labels = map_labs,
      limits = c(col_min, col_max), na.value = '#AAAAAA'
    ) +
    coord_sf(crs = india_refsys) +
    labs(
      title = titles[[ag]], x='', y='',
      fill=paste0(legend_ag,'\nMortality\nRate')) +
    theme_minimal() +
    theme(
      axis.text.x=element_blank(), axis.ticks.x=element_blank(),
      axis.text.y=element_blank(), axis.ticks.y=element_blank(),
      panel.grid.major = element_line(colour = 'transparent')
    )
  return(fig)
}

fig1_admin_maps <- list(
  make_admin_map('neonatal') + theme(plot.margin=unit(c(.25, 0, -.25, 0), "cm")),
  make_admin_map('infant') + theme(plot.margin=unit(c(0, 0, 0, 0), "cm")),
  make_admin_map('under5') + theme(plot.margin=unit(c(-.25, 0, .25, 0), "cm"))
)
fig1_layout_matrix <- matrix(1:3, ncol=1)
fig1_grobs <- lapply(fig1_admin_maps, ggplotGrob)

png(file.path(out_dir, 'fig1_mort_2000_2017.png'), height=9, width=7.5, units='in', res=map_res)
grid.arrange(grobs=fig1_grobs, layout_matrix=fig1_layout_matrix)
dev.off()



## ** (FIG 2: AROC by district in U5MR, 2000 to 2017) **

# Data prep
fig2_data <- copy(aroc$under5)
fig2_breaks <- c(-0.05, -0.025, -0.01, 0)
fig2_labs <- c('5%', '2.5%', '1%', 'No decline')
fig2_cols <- c("#400080", "#bd8cd9", "#FFFFE6", "#d9b38c", "#da730b")
fig2_data[ mean < min(fig2_breaks), mean := min(fig2_breaks) + 1E-5 ]
fig2_data[ mean > max(fig2_breaks), mean := max(fig2_breaks) - 1E-5 ]
fig2_sf <- merge(shps$ad2, fig2_data[, .(mean, ADM2_CODE)], by='ADM2_CODE', all.x=TRUE)

# Plotting
fig2 <- ggplot() +
  geom_sf(data = fig2_sf, aes(fill=mean), lwd=0.03, color=map_linecol) +
  geom_sf(data = shps$ad1, lwd=0.2, fill=NA, color=map_linecol) +
  scale_fill_gradientn(
    colors = fig2_cols, breaks = fig2_breaks, labels = fig2_labs,
    limits = range(fig2_breaks), na.value = '#AAAAAA'
  ) +
  coord_sf(crs = india_refsys) +
  labs(title = '', x='', y='', fill='Annualized\nRate of\nDecline,\n2000-2017') +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent')
  )

png(file.path(out_dir, 'fig2_u5m_aroc_2000_2017.png'), height=5, width=5, units='in', res=map_res)
print(fig2)
dev.off()



## ** (FIG 3: High-low plot of district U5MR, IMR, and NMR by Indian state) **

hilo_plot_year_start <- 2000
hilo_plot_year_end <- 2017

# Prepare data
hilo_data <- rbindlist(lapply(age_groups,
  function(x) (summ[[x]]
    [, age_group := titles[x]]
    [year %in% c(hilo_plot_year_start, hilo_plot_year_end)]
    [, mean_q := mean_q * 1000]
  )
))
hilo_data$ag_factor <- factor(hilo_data$age_group, levels=unlist(titles))
hilo_ad1 <- hilo_data[admin_level==1,]
hilo_ad2 <- hilo_data[admin_level==2,]
# Merge on admin1 code
mts <- c('Andaman & Nicobar Island', 'Chandigarh', 'Dadra & Nagar Haveli', 'Daman & Diu',
         'Lakshadweep', 'Puducherry')
hilo_ad2[ADM1_NAME %in% mts, ADM1_NAME := "The Six Minor Territories" ]
hilo_ad2[ADM1_NAME == "Jammu & Kashmir", ADM1_NAME := "Jammu and Kashmir" ]
hilo_ad2[ADM1_NAME == "Delhi", ADM1_NAME := "NCT of Delhi" ]
hilo_ad2[hilo_ad1[year==hilo_plot_year_start,], ADM1_CODE := i.ADM1_CODE, on='ADM1_NAME']
# Sort by mortality in final plotting year
hilo_plot_rank <- (hilo_ad1
  [year==hilo_plot_year_end & age_group =='Neonatal (<1 month)', ]
  [order(-mean_q)]
  [, plot_order := .I ]
)
hilo_ad1[hilo_plot_rank, plot_order := i.plot_order, on ='ADM1_CODE']
hilo_ad2[hilo_plot_rank, plot_order := i.plot_order, on ='ADM1_CODE']

years_palette <- c('#888888', 'blue')
names(years_palette) <- as.character(c(hilo_plot_year_start, hilo_plot_year_end))
hilo_ad2[, year_color := as.character(year) ]

# Prep figure
targets_dt <- data.table(age_group = titles, q_target = c(16, 28, 23))
targets_dt$ag_factor <- factor(targets_dt$age_group, levels=unlist(titles))

p_ly <- position_nudge(x=.2)
fig3 <- ggplot(
    data = hilo_ad1[year==hilo_plot_year_start,],
    aes(x=reorder(ADM1_NAME, plot_order), y=mean_q)
  ) +
  facet_wrap('ag_factor', ncol=1, scales='free_y') +
  # Point showing goal
  geom_hline(data = targets_dt, aes(yintercept = q_target), linetype=2, color='#888888') +
  # All admin2s
  geom_point(
    data=hilo_ad2[year==hilo_plot_year_start,], aes(color=year_color),
    size=1.1, shape=18, alpha=.5
  ) +
  geom_point(
    data=hilo_ad2[year==hilo_plot_year_end,], aes(color=year_color),
    size=1.1, shape=18, alpha=.5, position=p_ly
  ) +
  # Admin1 midpoints
  geom_point(data=hilo_ad1[year==hilo_plot_year_start,], shape=18, alpha=.5, size=2.2) +
  geom_point(data=hilo_ad1[year==hilo_plot_year_end,], shape=18, alpha=.5, size=2.2, position=p_ly) +
  scale_color_manual(values = years_palette) +
  guides(color = guide_legend(
    override.aes = list(size = 3), nrow = 1, title.position = "left"
  )) +
  labs(x='State or UT', y = 'Deaths per 1,000 live births', color='Year') +
  theme_bw() +
  theme(
    text = element_text(size=12),
    axis.text.x = element_text(size=9, angle = 60, hjust = 1),
    legend.position = c(.87, .96),
    legend.text = element_text(margin = margin(l = -8))
  )

png(glue('{out_dir}/fig3_hilo_all_ages.png'), height=8, width=7.5, units='in', res=map_res)
print(fig3)
dev.off()



## ** (FIG 4: U5MR in 2025 (projected) nationwide) **

# Data prep
fig4_data <- copy(proj$under5[year==2025, ])

target_2025 <- 0.023
fig4_labs <- breaks$under5
fig4_breaks <- fig4_labs / 1000

fig4_data[, target := 'Uncertain']
fig4_data[upper < target_2025, target := 'Yes']
fig4_data[lower > target_2025, target := 'No']

fig4_sf <- merge(
  x = shps$ad2, y=fig4_data[, .(ADM2_CODE, mean, target)], by='ADM2_CODE', all.x=TRUE
)
fig4_sf[is.na(fig4_sf$target), 'target'] <- 'Uncertain'

# Plot overall rate
fig_4a <- ggplot() +
  geom_sf(data = fig4_sf, aes(fill=mean), lwd=0.03, color=map_linecol) +
  geom_sf(data = shps$ad1, lwd=0.2, fill=NA, color=map_linecol) +
  scale_fill_gradientn(
    colors = color_scheme, breaks = fig4_breaks, labels = fig4_labs,
    limits = range(fig4_breaks), na.value = '#AAAAAA', oob=scales::squish
  ) +
  coord_sf(crs = india_refsys) +
  labs(
    title = 'Under-5 Mortality Rate, 2025 (projected)', x='', y='',
    fill='Under-5\nMortality Rate'
  ) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent')
  )

colors_4b <- c("#33cccc", "#d1476a", "#ffffff")
names(colors_4b) <- c("Yes", "No", "Uncertain")
fig_4b <- ggplot() +
  geom_sf(data = fig4_sf, aes(fill=target), lwd=0.03, color=map_linecol) +
  geom_sf(data = shps$ad1, lwd=0.2, fill=NA, color=map_linecol) +
  scale_fill_manual(values = colors_4b, na.value = '#AAAAAA') +
  coord_sf(crs = india_refsys) +
  labs(title = 'Projected to meet NHP 2025 target?', x='', y='', fill='') +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent')
  )

fig4_maps <- list(ggplotGrob(fig_4a), ggplotGrob(fig_4b))
fig4_layout_matrix <- matrix(1:2, ncol=1)

png(file.path(out_dir, 'fig4_u5m_proj_2025.png'), height=7, width=4.5, units='in', res=map_res)
grid.arrange(grobs=fig4_maps, layout_matrix=fig4_layout_matrix)
dev.off()


## ** (FIG 5: IMR for all three, at the most detailed level available) **

# Prepare IMR by "natural division" in 2017
imr_dist <- summ$infant[(admin_level == 2) & (year == 2017), ]
imr_dist[, pop_est := mean_d / mean_q ]
imr_dist[, `:=` (mean_q=mean_q*1000, lower_q=lower_q*1000, upper_q=upper_q*1000)]
nd_grouping_cols <- c('ADM1_CODE','nd_code','nd_name','nd_idx')
imr_nd <- merge(
  x=nd_dist, y=imr_dist[, .(ADM2_CODE, mean_q, lower_q, upper_q, pop_est)],
  by.x='DIST_ID', by.y='ADM2_CODE'
)
imr_nd <- imr_nd[
  , .(mean_q = weighted.mean(mean_q, w=pop_est),
      lower_q = weighted.mean(lower_q, w=pop_est),
      upper_q=weighted.mean(upper_q, w=pop_est)),
  by=nd_grouping_cols
]
# Merge on data from the SRS 2017
imr_nd[srs_tabs, srs_imr := i.est_imr, on = c('ADM1_CODE','nd_code')]
imr_nd <- imr_nd[!is.na(srs_imr), ]
imr_nd <- merge(imr_nd, unique(nd_dist[, .(ADM1_CODE, ADM1_NAME)]))
# Get ratio of survey:SRS and check for significance
imr_nd[, srs_svy_ratio := srs_imr / mean_q ]
sig_ids <- c('Significantly higher', 'Significantly lower', 'Overlapping UIs')
imr_nd[, significance := sig_ids[3]]
imr_nd[srs_imr < lower_q, significance := sig_ids[2]]
imr_nd[srs_imr > upper_q, significance := sig_ids[1]]
fwrite(imr_nd, file=file.path(out_dir, 'plug_nd_imr_comparison.csv'))

# Merge on spatial data
imr_nd_spatial <- merge(
  x = nd_sf_sub,
  y = imr_nd[, .(nd_idx, srs_imr, srs_svy_ratio, significance)]
)


## Plot sub-maps
fig_5_base <- ggplot(data=imr_nd_spatial) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent')
  )
# Sub-map 1: plot estimated IMR from SRS data
map_breaks <- map_labs <- breaks$infant
col_max <- max(map_breaks)
col_min <- min(map_breaks)
imr_nd_spatial[imr_nd_spatial$srs_imr > col_max ,'srs_imr'] <- col_max - 1E-5
imr_nd_spatial[imr_nd_spatial$srs_imr < col_min ,'srs_imr'] <- col_min + 1E-5

fig_5a <- fig_5_base +
  geom_sf(data = shps$ad1, lwd=0, fill='#AAAAAA', color=NA) +
  geom_sf(data = imr_nd_spatial, aes(fill=srs_imr), lwd=0.05, color='#444444', linetype=3) +
  geom_sf(data = shps$ad1, lwd=.3, fill=NA, color=map_linecol) +
  coord_sf(crs = india_refsys) +
  labs(
    title = 'A', x='', y='', fill='Infant\nMortality Rate\n(SRS)'
  ) +
  scale_fill_gradientn(
    colors = color_scheme, breaks = map_breaks, labels = map_labs,
    limits = c(col_min, col_max), na.value = '#AAAAAA'
  ) +
  theme(plot.margin=unit(c(.25, 0, -.75, 0), "cm"))

# Sub-map 2: plot ratio of SRS / survey data
fig_5b_colors <- RColorBrewer::brewer.pal(name='PiYG', n=11)
fig_5b_breaks <- c(.5, .75, 1., 1.25, 1.5)
fig_5b_labs <- c('0.5 (survey higher)', '0.75', '1.0', '1.25', '1.5 (SRS higher)')

fig_5b <- fig_5_base +
  geom_sf(data = shps$ad1, lwd=0, fill='#AAAAAA', color=NA) +
  geom_sf(data = imr_nd_spatial, aes(fill=srs_svy_ratio), lwd=0.05, color='#444444', linetype=3) +
  geom_sf(data = shps$ad1, lwd=.3, fill=NA, color=map_linecol) +
  coord_sf(crs = india_refsys) +
  labs(
    title = 'B', x='', y='', fill='Ratio of IMR\nestimates:\nSRS vs. Survey'
  ) +
  scale_fill_gradientn(
    colors = fig_5b_colors, breaks = fig_5b_breaks, labels = fig_5b_labs,
    limits = c(0.45, 1.55), na.value = '#AAAAAA'
  ) +
  theme(plot.margin=unit(c(-.25, 0, -.25, 0), "cm"))

# Sub-map 3: plot whether the difference is significant
fig_5c_colors <- c(fig_5b_colors[c(10, 2)], '#FFFFFF')
names(fig_5c_colors) <- sig_ids

fig_5c <- fig_5_base +
  geom_sf(data = shps$ad1, lwd=0, fill='#AAAAAA', color=NA) +
  geom_sf(data = imr_nd_spatial, aes(fill=significance), lwd=0.05, color='#444444', linetype=3) +
  geom_sf(data = shps$ad1, lwd=.3, fill=NA, color=map_linecol) +
  coord_sf(crs = india_refsys) +
  labs(
    title = 'C', x='', y='', fill='SRS estimate\nsignificantly\ndifferent?'
  ) +
  scale_fill_manual(
    values = fig_5c_colors,
    guide = guide_legend(override.aes = list(linetype=1))
  ) +
  theme(plot.margin=unit(c(-.75, 0, .25, 0), "cm"))

# Figure 5: Combine 'em
layout_matrix <- matrix(1:3, ncol=1)
fig5_grobs <- lapply(list(fig_5a, fig_5b, fig_5c), ggplotGrob)

png(file.path(out_dir, 'fig5_srs_comparison.png'), height=9, width=5.5, units='in', res=map_res)
grid.arrange(grobs = fig5_grobs, layout_matrix = layout_matrix)
dev.off()


## ** (FIG 6: Number of U5 individuals covered by year from survey data, 2000 to 2017) **

raw_data <- fread('C:/Users/nathenry/Desktop/india_sample_size.csv')
# Drop sub-year age bins
raw_data <- raw_data[!(ab %in% c('PNN1', 'PNN2')) & (country=='IND') & (year >= 2000), ]
sample_sizes <- rbind(
  raw_data[, .(N = sum(N*weight)), by=year][order(year)],
  data.table(year=2017, N=0)
)

fig6 <- ggplot(data=sample_sizes, aes(x=year, y=N)) +
  geom_line() +
  scale_x_continuous(limits=c(2000, 2017), breaks=c(seq(2000, 2015, by=5), 2017)) +
  scale_y_continuous(
    limits=c(0, 1.15E6), breaks=seq(0, 1E6, by=2.5E5), labels=c('0','250','500','750','1,000')
  ) +
  labs(
    title='Sample size for retrospective birth histories over time',
    x='Year', y='Sample size for children under 5\n(in thousands of person-years)'
  ) +
  theme_bw()

png(file.path(out_dir,'fig6_sample_size.png'), height=4, width=6, units='in', res=map_res)
print(fig6)
dev.off()


## ** (FIG 7: Ranked districts in 2000, 2017, and 2025 (est) with uncertainty showing
##      flattened outcomes in 2017) **

u5_past <- summ$under5[(admin_level == 2) & (year %in% c(2000, 2017)), ]
u5_future <- proj$under5[year==2025, ]
qstats <- c('mean','lower','upper')
setnames(u5_future, qstats, paste0(qstats,'_q'))

keep_fields <- c('ADM1_NAME','ADM2_NAME','ADM2_CODE','year','mean_q','upper_q','lower_q')
f7_data <- rbindlist(
  list(u5_past[, ..keep_fields], u5_future[, ..keep_fields]), use.names=T, fill=T
)[order(ADM1_NAME, ADM2_NAME, year)]
f7_data[, `:=` (mean_q = mean_q*1E3, lower_q=lower_q*1E3, upper_q=upper_q*1E3)]

# Name updates
f7_data[ADM2_NAME == 'Chittaurgarh', ADM2_NAME := 'Chittorgarh']
f7_data[ADM2_NAME == 'Dhaulpur', ADM2_NAME := 'Dholpur']
f7_data[ADM2_NAME == 'Jalor', ADM2_NAME := 'Jalore']
f7_data[ADM2_NAME == 'Jhunjhunun', ADM2_NAME := 'Jhunjhunu']
f7_data[ADM2_NAME == 'Sri Potti Sriramulu Nellore', ADM2_NAME := 'Nellore']
f7_data[ADM2_NAME == 'Y.S.R', ADM2_NAME := 'Y.S.R.']
f7_data[ADM2_NAME == 'Vizianagaram', ADM2_NAME := 'Vizianagram']
f7_data[ADM2_NAME == 'Mahbubnagar', ADM2_NAME := 'Mahabubnagar']

# Rank them
f7_ranks <- f7_data[year==2000,][, plot_rank := frank(-mean_q), by=ADM1_NAME]
f7_data[f7_ranks, plot_rank := i.plot_rank, on=c('ADM1_NAME','ADM2_NAME')]
f7_data[, year := as.character(year)]

# Make plots
f7_sub_dt <- f7_data[ADM1_NAME %in% c('Andhra Pradesh','Rajasthan','Telangana'), ]
fig_7 <- ggplot(data = f7_sub_dt, aes(x=reorder(ADM2_NAME, plot_rank), y=mean_q, color=year)) +
  facet_wrap('ADM1_NAME', ncol=1, scales='free_x') +
  geom_crossbar(
    aes(ymin=lower_q, ymax=upper_q), position=position_dodge(width=.5),
    width=0, show.legend = FALSE) +
  geom_point(shape=18, size=1.5, position=position_dodge(width=.5)) +
  labs(title=NULL, x='District (Ranked)', y='Under-5 Mortality Rate', color='Year') +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  theme_bw() +
  theme(
    text = element_text(size=16),
    axis.text.x = element_text(size=7, angle = 60, hjust = 1),
    legend.text = element_text(margin = margin(l = -8))
  ) +
  geom_hline(yintercept = 25, linetype = 2, color = '#888888')

png(file.path(out_dir, 'fig7_rankings.png'), height=9, width=7, units='in', res=map_res*1.5)
print(fig_7)
dev.off()
