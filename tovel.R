# Project: Tovel Lake Modeling Study
# Description: General Lake Model (GLM) for Lake Tovel
# Last modified: 27 December 2021 
# Author: Robert Ladwig, Ulrike Obertegger

#### Setup ####
cat("\f")
rm(list = ls())

# if you're using Rstudio:
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# we will need these packages
library(glmtools)
library(GLM3r)
library(rLakeAnalyzer)
library(tidyverse)
library(lubridate)
library(gotmtools)

library(ncdf4)
# Now open the file and read its data
nc <- nc_open('output/output.nc')
data <- ncvar_get(nc)
print("here is the data in the file:")
print(data)
nc_close( nc )


hypso <- read_csv('bathymetry.csv')

# check out which R version we're currently using
glm_version() #version 3.0.5

#### Example 1: reading the namelist file into R  ####
# glm_template = 'glm3-template.nml' 
sim_folder <- getwd()
out_file <- file.path(sim_folder, "output","output.nc")
# field_data <- file.path(sim_folder,"LakeEnsamble_WaterTemperature_2020.csv")
field_data <- file.path(sim_folder,"bc_feb22/LakeEnsamble_WaterTemperature_2017-2020.csv")
nml_file <- file.path(sim_folder, 'glm3.nml')

# read initial field data
field <- read_csv(field_data) %>%
  filter(Time == min(Time))

# read example configuration into memory
eg_nml <- read_nml(nml_file = file.path(sim_folder,'glm3.nml'))

eg_nml <- set_nml(eg_nml, 'H', hypso$Depth_meter - min(hypso$Depth_meter) + (1178 - max(hypso$Depth_meter)))
eg_nml <- set_nml(eg_nml, 'A', rev(hypso$Area_meterSquared))
eg_nml <- set_nml(eg_nml, 'bsn_vals', length(hypso$Area_meterSquared))
eg_nml <- set_nml(eg_nml, 'crest_elev', 1178)
eg_nml <- set_nml(eg_nml, 'num_depths', length(field$Depth_meter))
eg_nml <- set_nml(eg_nml, 'the_depths', field$Depth_meter)
eg_nml <- set_nml(eg_nml, 'the_temps', field$Water_Temperature_celsius)
eg_nml <- set_nml(eg_nml, 'the_sals', rep(0, length(field$Depth_meter)))

# remotes::install_github("aemon-j/gotmtools", force = TRUE)
# meteo <- read_csv('bc_meteo.csv')
meteo <-c()
# orig_meteo <- read_csv('LakeEnsamble_Meteo_2020.csv')
orig_meteo <- read_csv('bc_feb22/LakeEnsamble_Meteo_2017-2020.csv')



cc <- gotmtools::calc_cc(date = orig_meteo$datetime,
                   airt = orig_meteo$Air_Temperature_celsius,
                   relh = orig_meteo$Relative_Humidity_percent,
                   swr = orig_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared,
                   lat = eg_nml$morphometry$latitude, lon = eg_nml$morphometry$longitude,
                   elev = eg_nml$morphometry$crest_elev)

lw <- gotmtools::calc_in_lwr(cc = cc,
                       airt = orig_meteo$Air_Temperature_celsius,
                       relh = orig_meteo$Relative_Humidity_percent)
meteo$Date <- orig_meteo$datetime
meteo$ShortWave <- orig_meteo$Shortwave_Radiation_Downwelling_wattPerMeterSquared
meteo$LongWave <- lw
meteo$Clouds <- cc
meteo$WindSpeed <- orig_meteo$Ten_Meter_Elevation_Wind_Speed_meterPerSecond
meteo$AirTemp <- orig_meteo$Air_Temperature_celsius
meteo$RelHum <- orig_meteo$Relative_Humidity_percent
meteo$Rain <- ifelse(meteo$AirTemp > 0, orig_meteo$Precipitation_millimeterPerDay/1000, 0)
meteo$Snow <- ifelse(meteo$AirTemp <= 0, orig_meteo$Precipitation_millimeterPerDay/1000, 0)
meteo <- as.data.frame(meteo)

out <- boxplot.stats(meteo$AirTemp)$out
out_ind <- which(meteo$AirTemp %in% c(out))
meteo$AirTemp[out_ind] = NA
meteo$AirTemp <- zoo::na.approx(meteo$AirTemp)

write_csv(meteo, 'bc_meteo.csv')

m.meteo <- reshape2::melt(meteo, id = 'Date')
ggplot(m.meteo) +
  geom_line(aes(Date, value)) +
  facet_wrap(~ variable, scales = 'free') +
  theme_minimal()

obs_wlvl <- read_csv('bc_feb22/waterlevel_2017-2020.csv')
eg_nml <- set_nml(eg_nml, 'lake_depth', obs_wlvl$`waterlevel [m]`[1] - (1178 - max(hypso$Depth_meter)))

if (max(get_nml_value(eg_nml, 'the_depths')) > get_nml_value(eg_nml, 'lake_depth')){
  dat = field$Depth_meter
  idx <- which.max(dat)
  dat[20] = floor( get_nml_value(eg_nml, 'lake_depth'))
  eg_nml <- set_nml(eg_nml, 'the_depths', dat)
} 

# write modified values into namelist file
write_nml(eg_nml, file = nml_file)


# light
light <- c()
df <- read_csv('bc_feb22/Secchi_Tovel_2017-2020.csv')
light$date <- as.POSIXct(paste0(df$year,'-',formatC(df$month, flag=0, width =2),'-',df$day,' 00:00:00'))
light$Kd <-2/ df$`secchi disk [m]`

light <- as.data.frame(light)

write_csv(light, 'bc_light.csv')


q <- rep(0, nrow(meteo))
h <- rep(0, nrow(meteo))
h[1] <- obs_wlvl$waterlevel[1]
for (i in 2:length(q)){
  q[i] <- ((obs_wlvl$waterlevel[i] - obs_wlvl$waterlevel[i-1])/3600 - meteo$Rain[i]/86400 - meteo$Snow[i]/86400) * max(hypso$Area_meterSquared)+ 0.2
  if (meteo$AirTemp[i] <= 0 & q[i] > 0 & month(meteo$Date[i] ) < month(4)){
    q[i] = 0
  }
  if (meteo$AirTemp[i] <= 0 & q[i] < 0 & month(meteo$Date[i] ) < month(4)){
    q[i] = q[i]/2
  }
  h[i] <-  h[i-1] + 3600 * (q[i]/max(hypso$Area_meterSquared) + meteo$Rain[i]/86400 + meteo$Snow[i]/86400 - 0.2/max(hypso$Area_meterSquared))
}
idx <- which(q > 0)
idy <- which(q < 0)
q_in = q
q_out = q
q_in[idy] = 0 
q_out[idx] = 0


# write_csv(data.frame('time' = meteo$Date, 'FLOW' = (q_in), 'SALT' = rep(0, length(q)), 'TEMP' = rep(5, length(q))), 'bc_inflow.csv')

# write_csv(data.frame('time' = meteo$Date, 'FLOW' = (q_in), 'SALT' = rep(0, length(q)), 'TEMP' = meteo$AirTemp), 'bc_inflow.csv')
# write_csv(data.frame('time' = meteo$Date, 'FLOW' = abs(q_out)), 'bc_outflow.csv')

ggplot(data.frame('time' = meteo$Date, 'FLOW' = (q_in), 'SALT' = rep(0, length(q)), 'TEMP' = meteo$AirTemp))+
  geom_line(aes(time, FLOW, col = 'FLOW'))+
  geom_point(aes(time, TEMP, col = 'TEMP'), alpha = 0.2) +
  theme_bw()

inf <- read_csv('bc_feb22/Q_Pozzol_2017-2020.csv')
inf$datetime <- as.POSIXct(paste0(inf$date,' 12:00:00'))
meteo <- read_csv('bc_meteo.csv')
meteo$Date <- as.POSIXct(meteo$Date)
inf$temp <- meteo$AirTemp[match(inf$datetime, meteo$Date)]
df <- data.frame(time = inf$datetime,
                 FLOW = inf$`inflow Q [m3/s]`,
                 SALT = 0,
                 TEMP = inf$temp)
write_csv(df, 'bc_inflow.csv')

plot(q_in)
# seepage of 200 l/s or 0.2 m3/s
0.2 / max(hypso$Area_meterSquared)

# run GLM
GLM3r::run_glm(sim_folder, verbose = T)

field_data <- read_csv("bc_feb22/LakeEnsamble_WaterTemperature_2017-2020.csv")
split_hour = as.POSIXct(cut(as.POSIXct(field_data$Time), breaks = "60 mins"), format = '%Y-%m-%d %H:%M:%S', tz = 'UTC') # summrise given mins
field_data$hour = split_hour # make hourly vaiable
ag = aggregate(. ~ hour, field_data %>% group_by(Depth_meter, hour), mean)

ag = field_data[which(field_data$Time %in% split_hour), ]

write_csv(data.frame('DateTime' = as.POSIXct(ag$Time, format = '%Y-%m-%d %H:%M:%S', tz = 'UTC'),'Depth' = ag$Depth_meter, 'temp' = ag$Water_Temperature_celsius
), 'bc_field.csv', )

# calibrate Tovel

glm_template = 'glm3-template.nml'
file.copy(glm_template, 'glm3.nml', overwrite = TRUE)
file.remove('glm4.nml')

var = 'temp'         # variable to which we apply the calibration procedure
path = getwd()       # simulation path/folder
nml_file = nml_file  # path of the nml configuration file that you want to calibrate on
glm_file = nml_file # # path of the gml configuration file
# which parameter do you want to calibrate? a sensitivity analysis helps
calib_setup <- data.frame('pars' = as.character(c('wind_factor','lw_factor','ch',
                                                  'coef_mix_hyp', 'coef_mix_KH',
                                                  'sw_factor','cd')),
                          'lb' = c(0.7,0.7,0.0007,1e-20, 0.1,0.7,5e-4),
                          'ub' = c(1.3,1.3,0.002,1e-5, 0.9,1.3,0.002),
                          'x0' = c(1,1,0.0013,1e-10, 0.3,1,0.0013))
print(calib_setup)

# calib_setup <- data.frame('pars' = as.character(c('wind_factor','lw_factor','ch',
#                                                   'coef_mix_hyp', 'coef_mix_KH')),
#                           'lb' = c(0.7,0.7,5e-4,1e-8, 0.1),
#                           'ub' = c(1.3,1.3,0.002,1e-1, 0.9),
#                           'x0' = c(1,1,0.0013,1e-5, 0.5))
# print(calib_setup)
glmcmd = NULL        # command to be used, default applies the GLM3r function
# glmcmd = '/Users/robertladwig/Documents/AquaticEcoDynamics_gfort/GLM/glm'        # custom path to executable
# Optional variables
first.attempt = TRUE # if TRUE, deletes all local csv-files that stores the 
#outcome of previous calibration runs
period = get_calib_periods(nml_file, ratio = 3) # define a period for the calibration, 
# this supports a split-sample calibration (e.g. calibration and validation period)
# the ratio value is the ratio of calibration period to validation period
print(period)

period$total$start <- as.character('2017-01-01 00:00:00')
period$calibration$start <- as.character('2017-01-01 00:00:00')
period$calibration$stop <- as.character('2020-12-01 23:00:00')
period$validation$start <- as.character('2020-11-01 00:00:00')
period$validation$stop <- as.character('2020-12-01 23:00:00')
str(period)

scaling = TRUE       # scaling of the variables in a space of [0,10]; TRUE for CMA-ES
verbose = TRUE
method = 'CMA-ES'    # optimization method, choose either `CMA-ES` or `Nelder-Mead`
metric = 'RMSE'      # objective function to be minimized, here the root-mean square error
target.fit = 1.0     # refers to a target fit of 2.0 degrees Celsius (stops when RMSE is below that)
target.iter = 300    # refers to a maximum run of 20 calibration iterations (stops after that many runs)
plotting = TRUE      # if TRUE, script will automatically save the contour plots
output = out_file    # path of the output file
field_file = 'bc_feb22/LakeEnsamble_WaterTemperature_2017-2020_update.csv'#'bc_field.csv' # path of the field data
conversion.factor = 1 # conversion factor for the output, e.g. 1 for water temp.

calibrate_sim_2d(var = 'temp', path = getwd(), 
              field_file = field_file, 
              nml_file = nml_file, 
              glm_file = glm_file, 
              calib_setup = calib_setup, 
              glmcmd = NULL, first.attempt = TRUE, 
              period = period, 
              scaling = TRUE, method = 'CMA-ES', metric = 'RMSE', 
              target.fit = target.fit, target.iter = target.iter, 
              plotting = TRUE, 
              output = output, 
              verbose = TRUE,
              conversion.factor = 1)

# 1D
glm_template = 'glm3-template.nml'
file.copy(glm_template, 'glm3.nml', overwrite = TRUE)
file.remove('glm4.nml')

var = 'temp'         # variable to which we apply the calibration procedure
path = getwd()       # simulation path/folder
nml_file = nml_file  # path of the nml configuration file that you want to calibrate on
glm_file = nml_file # # path of the gml configuration file
# which parameter do you want to calibrate? a sensitivity analysis helps
calib_setup <- data.frame('pars' = as.character(c('outflow_factor', 'seepage_rate','inflow_factor', 'rain_factor')),
                          'lb' = c(0.1, -0.12, 0.1,0.7), 
                          'ub' = c(2, 0.01, 2,1.3),
                          'x0' = c(1, -0.005, 1,1))
print(calib_setup)

# calib_setup <- data.frame('pars' = as.character(c('wind_factor','lw_factor','ch',
#                                                   'coef_mix_hyp', 'coef_mix_KH')),
#                           'lb' = c(0.7,0.7,5e-4,1e-8, 0.1),
#                           'ub' = c(1.3,1.3,0.002,1e-1, 0.9),
#                           'x0' = c(1,1,0.0013,1e-5, 0.5))
# print(calib_setup)
glmcmd = NULL        # command to be used, default applies the GLM3r function
# glmcmd = '/Users/robertladwig/Documents/AquaticEcoDynamics_gfort/GLM/glm'        # custom path to executable
# Optional variables
first.attempt = TRUE # if TRUE, deletes all local csv-files that stores the 
#outcome of previous calibration runs
period = get_calib_periods(nml_file, ratio = 3) # define a period for the calibration, 
# this supports a split-sample calibration (e.g. calibration and validation period)
# the ratio value is the ratio of calibration period to validation period
print(period)

period$total$start <- as.character('2017-01-01 00:00:00')
period$calibration$start <- as.character('2017-01-01 00:00:00')
period$calibration$stop <- as.character('2020-12-01 23:00:00')
period$validation$start <- as.character('2020-11-01 00:00:00')
period$validation$stop <- as.character('2020-12-01 23:00:00')
str(period)

scaling = TRUE       # scaling of the variables in a space of [0,10]; TRUE for CMA-ES
verbose = TRUE
method = 'CMA-ES'    # optimization method, choose either `CMA-ES` or `Nelder-Mead`
metric = 'RMSE'      # objective function to be minimized, here the root-mean square error
target.fit = 0.1   # refers to a target fit of 2.0 degrees Celsius (stops when RMSE is below that)
target.iter = 1000    # refers to a maximum run of 20 calibration iterations (stops after that many runs)
plotting = TRUE      # if TRUE, script will automatically save the contour plots
output = out_file    # path of the output file
field_file = 'bc_feb22/wlvl.csv' # path of the field data
conversion.factor = 1 # conversion factor for the output, e.g. 1 for water temp.

calibrate_sim_1d(var = 'waterlevel', path = getwd(), 
              field_file = field_file, 
              nml_file = nml_file, 
              glm_file = glm_file, 
              calib_setup = calib_setup, 
              glmcmd = NULL, first.attempt = TRUE, 
              period = period, 
              scaling = TRUE, method = 'CMA-ES', metric = 'RMSE', 
              target.fit = target.fit, target.iter = target.iter, 
              plotting = TRUE, 
              output = output, 
              verbose = TRUE,
              conversion.factor = 1)


# run GLM
GLM3r::run_glm(sim_folder, verbose = T)
g <- plot_var_compare(nc_file = out_file,
                 field_file = 'bc_field.csv',zlim = c(-1,35),
                 var_name = 'temp', precision = 'days', obs.alpha = 0.1,
                 obs.size = 0.1); g
ggsave(plot = g,'heatmap.png', dpi = 300, width = 15, height = 10, units = 'in')

# visualize change of water table over time
# obs_wlvl <- read_csv('bc_feb22/waterlevel_2017-2020.csv')
obs_wlvl <- read_csv('bc_feb22/wlvl.csv')
water_height <- get_surface_height(file = out_file) 
# water_height$surface_height = water_height$surface_height +  (1178 - max(hypso$Depth_meter))
g.wlvl <- ggplot(water_height, aes(DateTime, surface_height , col = 'sim')) +
  geom_line() +
  geom_point(data = obs_wlvl, aes(DateTime,  `waterlevel [m]`, col = 'obs')) +
  # geom_line(data = data.frame('time' = meteo$Date, 'h' = h), aes(time, h, col = 'approx'))+
  ggtitle('Surface water level') +
  # ylim(1170,1180) +
  xlab(label = '') + ylab(label = 'Water level (m)') +
  theme_bw();g.wlvl
ggsave('waterlvl.png', dpi = 300, width = 15, height = 10, units = 'in')

idx <- match(water_height$DateTime, obs_wlvl$DateTime)
res <- water_height$surface_height - obs_wlvl$`waterlevel [m]`[idx]
print(sqrt(sum(res^2, na.rm = T)/length(obs_wlvl$`waterlevel [m]`[idx])))

ggplot(data.frame('time' = obs_wlvl$datetime, 'sim' = c(1177,water_height$surface_height),
                  'obs' = obs_wlvl$`waterlevel [m]`), aes(obs, obs-sim)) +
  geom_point() +
  ggtitle('Surface water level') +
  xlab(label = '') + ylab(label = 'Water level (m)') +
  theme_minimal()

# visualize ice formation over time
ice_thickness <- get_ice(file = out_file)
g.ice <- ggplot(ice_thickness, aes(DateTime, `ice(m)`)) +
  geom_line() +
  ggtitle('Ice') +
  xlab(label = '') + ylab(label = 'Ice thickness (m)') +
  theme_bw();g.ice
ggsave('ice_thick.png', dpi = 300, width = 15, height = 10, units = 'in')

# visualize change of surface water temp. over time
surface_temp <- get_var(file = out_file, 
                        var_name = 'temp',
                        reference = 'surface',
                        z_out = 1)

middle_temp <- get_var(file = out_file, 
                        var_name = 'temp',
                        reference = 'surface',
                        z_out = 15)

# visualize change of bottom water temp. over time
bottom_temp <- get_var(file = out_file, 
                       var_name = 'temp',
                       reference = 'surface',
                       z_out = 30)

test = read_csv(field_data) 
obs_surf <- read_csv(file.path(sim_folder,"bc_feb22/LakeEnsamble_WaterTemperature_2017-2020_update.csv")) %>%
  filter(Depth == 1)#
obs_middle <- read_csv(file.path(sim_folder,"bc_feb22/LakeEnsamble_WaterTemperature_2017-2020_update.csv")) %>%
  filter(Depth == 15)#
obs_surftest <- read_csv('bc_field.csv') %>%
  filter(Depth == 1)
obs_bottom <- read_csv(file.path(sim_folder,"bc_feb22/LakeEnsamble_WaterTemperature_2017-2020_update.csv")) %>%
  filter(Depth == 30)

ggplot(surface_temp, aes(DateTime, temp_1)) +
  geom_point(data = obs_surf, aes(DateTime, temp, col = 'obs surface'), alpha = 0.2) +
  geom_point(data = obs_middle, aes(DateTime, temp,  col = 'obs middle'), alpha = 0.2) +
  geom_point(data = obs_bottom, aes(DateTime, temp,  col = 'obs bottom'), alpha = 0.2) +
  geom_line(aes(col = 'sim surface')) +
  geom_line(data = middle_temp, aes(DateTime, temp_15, col = 'sim middle')) +
  geom_line(data = bottom_temp, aes(DateTime, temp_30, col = 'sim bottom')) +
  ggtitle('Water temperature') +
  xlab(label = '') + ylab(label = 'Temp. (deg C)') +
  theme_bw()
ggsave('watertemp.png', dpi = 300, width = 15, height = 10, units = 'in')

g.surf <- ggplot(surface_temp, aes(DateTime, temp_1)) +
  geom_point(data = obs_surf, aes(DateTime, temp,  col = 'obs surface'), alpha = 0.2) +
  # geom_point(data = obs_middle, aes(Time, Water_Temperature_celsius, col = 'obs middle'), alpha = 0.2) +
  # geom_point(data = obs_bottom, aes(Time, Water_Temperature_celsius, col = 'obs bottom'), alpha = 0.2) +
  geom_line(aes(col = 'sim surface')) +
  # geom_line(data = middle_temp, aes(DateTime, temp_15, col = 'sim middle')) +
  # geom_line(data = bottom_temp, aes(DateTime, temp_30, col = 'sim bottom')) +
  ggtitle('Surface water temperature: 1 m') +
  xlab(label = '') + ylab(label = 'Temp. (deg C)') +
  theme_bw();g.surf


g.mid <- ggplot(surface_temp, aes(DateTime, temp_1)) +
  # geom_point(data = obs_surf, aes(Time, Water_Temperature_celsius, col = 'obs surface'), alpha = 0.2) +
  geom_point(data = obs_middle, aes(DateTime, temp,  col = 'obs middle'), alpha = 0.2) +
  # geom_point(data = obs_bottom, aes(Time, Water_Temperature_celsius, col = 'obs bottom'), alpha = 0.2) +
  # geom_line(aes(col = 'sim surface')) +
  geom_line(data = middle_temp, aes(DateTime, temp_15, col = 'sim middle')) +
  # geom_line(data = bottom_temp, aes(DateTime, temp_30, col = 'sim bottom')) +
  ggtitle('Middle water temperature: 15 m') +
  xlab(label = '') + ylab(label = 'Temp. (deg C)') +
  theme_bw();g.mid

g.bot <- ggplot(surface_temp, aes(DateTime, temp_1)) +
  # geom_point(data = obs_surf, aes(Time, Water_Temperature_celsius, col = 'obs surface'), alpha = 0.2) +
  # geom_point(data = obs_middle, aes(Time, Water_Temperature_celsius, col = 'obs middle'), alpha = 0.2) +
  geom_point(data = obs_bottom, aes(DateTime, temp, col = 'obs bottom'), alpha = 0.2) +
  # geom_line(aes(col = 'sim surface')) +
  # geom_line(data = middle_temp, aes(DateTime, temp_15, col = 'sim middle')) +
  geom_line(data = bottom_temp, aes(DateTime, temp_30, col = 'sim bottom')) +
  ggtitle('Bottom water temperature: 30 m') +
  xlab(label = '') + ylab(label = 'Temp. (deg C)') +
  theme_bw();g.bot

library(patchwork)
g.wlvl / g.ice / g.surf / g.mid / g.bot
ggsave('calibration.png', dpi = 300, width = 15, height = 10, units = 'in')



temp_rmse <- compare_to_field(nc_file = out_file, 
                              field_file = 'bc_field.csv',
                              metric = 'water.temperature', 
                              as_value = FALSE, 
                              precision= 'hours')
print(temp_rmse)

plot_var(nc_file = out_file, 
         var_name = 'temp',
         cex=3)
plot_var_compare(nc_file = out_file,
                 field_file = 'bc_field.csv',
                 var_name = 'temp', precision = 'days', obs.alpha = 0.1,
                 obs.size = 0.1)
ggsave('contour.png', dpi = 300, width = 15, height = 10, units = 'in')



# lakefile <- read_csv('output/lake.csv')
# m.lakefile <- reshape2::melt(lakefile, id = 'time')
# ggplot(m.lakefile) +
#   geom_line(aes(time, value))+
#   facet_wrap(~variable, scales = 'free') +
#   theme_minimal()
