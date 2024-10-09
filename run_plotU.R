here::i_am("./run_plotU.R")
library(renv)
library(here)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

Sys.setenv(TAR_PROJECT = "plotU")
tar_outdated(script = "_targets_plotU.R", store = "_targets_plotU")
# run in serial
system.time(tar_make(script = "_targets_plotU.R", store = "_targets_plotU"))
# run in parallel
system.time(tar_make_future(workers = 4L, script = "_targets_plotU.R", store = "_targets_plotU"))

# write graphics to file
for (region in v_region[1:4]) {
    target_name <- paste0("l_p_plotU_", region)
    l_p <- tar_read(eval(target_name), store = "_targets_plotU")
    fname <- here(paste0("output/output_", region), "r_uu_smp1_1000m.pdf")
    print(fname)
}

l_p <- tar_read(l_p_plotU_en, store = "_targets_plotU")
fname <- here("output/output_en", "r_uu_smp1_10000m.pdf")
l_p <- tar_read(l_p_plotU_sc, store = "_targets_plotU")
fname <- here("output/output_sc", "r_uu_smp1_10000m.pdf")
l_p <- tar_read(l_p_plotU_wa, store = "_targets_plotU")
fname <- here("output/output_wa", "r_uu_smp1_10000m.pdf")
l_p <- tar_read(l_p_plotU_ni, store = "_targets_plotU")
fname <- here("output/output_ni", "r_uu_smp1_10000m.pdf")

    pdf(fname)
        l_p$p_u_t
        l_p$p_u_d
        l_p$p_v1
        l_p$p_v2
    dev.off()

fname <- "/gws/nopw/j04/ceh_generic/plevy/luc_track/data-raw/pre1950/POPCEN/dt_A_DevA_km2.qs"
dt <- qread(fname)
dt <- drop_units(dt)
p <- ggplot(dt, aes(Year, urban))
p <- p + geom_line()
p <- p + facet_wrap(~ country, scale = "free_y")
p
dev.off()

r <- raster("/gws/nopw/j04/ceh_generic/plevy/luc_track/data-raw/LCM/Level1/r_U_lcm_100m_2019.tif")
r_en <- spCEH::maskByCountry(r, "England")
r_sc <- spCEH::maskByCountry(r, "Scotland")
r_wa <- spCEH::maskByCountry(r, "Wales")
r_ni <- spCEH::maskByCountry(r, "Northern Ireland")
cellStats(r == 5, sum)
cellStats(r_en == 5, sum)
cellStats(r_sc == 5, sum)
cellStats(r_wa == 5, sum)
cellStats(r_ni == 5, sum)

Year country woods crops grass rough     urban other
x <-
c( 2019,  "England",    NA,    NA,    NA,    NA, 10155,    NA,
   2019,    "Wales",    NA,    NA,    NA,    NA,   670,    NA,
   2019, "Scotland",    NA,    NA,    NA,    NA,  1270,    NA,
   2019,       "NI",    NA,    NA,    NA,    NA,   421,    NA)
x <- matrix(x, nrow = 4, ncol = 8, byrow = TRUE)
x <- as.data.table(matrix(x, nrow = 4, ncol = 8, byrow = TRUE))
names(x) <- names(dt)
dt <- rbind(dt, x)
dt[, urban := as.numeric(urban)]
dt[, Year := as.numeric(Year)]
quit(save = "no")
