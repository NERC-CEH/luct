# sbatch --account=short4hr slurm/run_sample_Upost_en.job
# sbatch --account=short4hr slurm/run_sample_Upost_sc.job
# sbatch --account=short4hr slurm/run_sample_Upost_wa.job
# sbatch --account=short4hr slurm/run_sample_Upost_ni.job
# sbatch --account=short4hr slurm/run_sample_Upost_uk.job
# R CMD BATCH --no-restore --no-save "--args  1 10000 en" slurm/sample_Upost.R "output/output_en/console_sample_Upost1.Rout"  > output/output_en/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  1 10000 sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost1.Rout"  > output/output_sc/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  1 10000 wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost1.Rout"  > output/output_wa/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  1 10000 ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost1.Rout"  > output/output_ni/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  1 1000  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost1.Rout"  > output/output_en/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  1 1000  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost1.Rout"  > output/output_sc/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  1 1000  wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost1.Rout"  > output/output_wa/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  1 1000  ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost1.Rout"  > output/output_ni/sample_Upost1.out &

# R CMD BATCH --no-restore --no-save "--args  1 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost1.Rout"  > output/output_en/sample_Upost1.out &
# R CMD BATCH --no-restore --no-save "--args  2 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost2.Rout"  > output/output_en/sample_Upost2.out &
# R CMD BATCH --no-restore --no-save "--args  3 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost3.Rout"  > output/output_en/sample_Upost3.out &
# R CMD BATCH --no-restore --no-save "--args  4 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost4.Rout"  > output/output_en/sample_Upost4.out &

# R CMD BATCH --no-restore --no-save "--args  1 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_1.Rout"  > output/output_sc/sample_Upost_100m_1.out &
# R CMD BATCH --no-restore --no-save "--args  2 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_2.Rout"  > output/output_sc/sample_Upost_100m_2.out &
# R CMD BATCH --no-restore --no-save "--args  3 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_3.Rout"  > output/output_sc/sample_Upost_100m_3.out &
# R CMD BATCH --no-restore --no-save "--args  4 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_4.Rout"  > output/output_sc/sample_Upost_100m_4.out &

#### pseudo-code
# read 4 inputs
# 1. current u and age dt_u_age.qs
# 2. static likelihood maps  L1-6
# 3. life tables a_iacs <- readRDS(here("data-raw/LifeTables", "IACS_LifeTables.rds"))
# 4. posterior Beta matrices a_B_post
# for time {
  # 1. read L_static
  # 2. update age_rev
  # 3. lookup life table
  # 4. update L1-6 = L_stat * L_dyn
  # for ij {
    # sample U_t
    # U_ch = append(U, U_t)
  # }
 # }

rm(list = ls(all = TRUE))
source("./_targets_packages.R")
library(targets)
library(here)
library(stringi)
library(wrswoR)
library(rasterVis)
library(ggthemes)
source("R/luc_track.R")
source("R/luct.R")
rasterOptions(tmpdir = "/work/scratch-nopw/bkruijt/raster")

set.seed(586)

# read arguments
# which sample in a_B_post to process
i_sample <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
if (is.na(i_sample)) i_sample <- 1
i_sample

# resolution to use in output raster (m)
res <- as.numeric(commandArgs(trailingOnly = TRUE))[2]
if (is.na(res)) res <- 10000
res

# region to restrict sampling to
region <- commandArgs(trailingOnly = TRUE)[3]
if (is.na(region)) region <- "uk"
region
i_region <- match(region, v_region) # 1:4
if (region == "uk") i_region <- 1:4
dir_output <- paste0("output/output_", region)

# define function to return prediction of v_U given lamda
getUt <- function(v_U_t, m_B, dt_L){
#v_U_t = dt[, u]
  diag(m_B) <- 0
  m_lamda <- as.matrix(dt_L)
# apply(dt_L, 2, summary)
# apply(m_lamda, 2, summary)
  #dim(m_lamda)
  #summary(m_lamda[,i_u])
  n_s <- length(v_U_t)
  # assume we go backwards, predicting t-1 from t
  v_U_tm1 <- v_U_t # initialise t-1 the same as t
# v_U_t  [v_ind]
# v_U_tm1[v_ind]
# v_mask [v_ind]
# m_lamda[v_ind, i_u]
# p = m_lamda[,i_u] * v_mask
# p[v_ind]
# summary(p)
# sum(p, na.rm = TRUE)
# min(m_lamda[,i_u], na.rm = TRUE)
  for (i_u in 1:n_u){
    for (j_u in 1:n_u){
      n_fromItoJ <- m_B[i_u, j_u] # check m_B defined correctly
      if (n_fromItoJ > 0){ 
        # mask: only those which are j_u at time t are possible candidates
        v_mask  <- v_U_t == j_u
        if (sum(v_mask, na.rm = TRUE) <= 0) next # no candidate cells
        v_mask[is.na(v_mask)]  <- 0
        # table(m_lamda[,i_u], v_mask)
        p <- m_lamda[,i_u] * v_mask
        # stopifnot is more rigorous, but seems to be single occurrence in Wales
        #stopifnot(sum(p, na.rm = TRUE) > 0)
        # replacing with a warning lets us carry on here
        if (sum(p, na.rm = TRUE) <= 0) warning(paste("No cells with positive probability going from", names_u[i_u], "to", names_u[j_u], "in", v_times[i_t]))                
        v_ind <- sample_int_crank(n_s, n_fromItoJ, prob = p)
        #m_lamda[v_ind,i_u]
        #sum(m_lamda[,i_u] > 0, na.rm = TRUE)
        #sort(v_ind)
        v_U_tm1[v_ind] <- i_u
        # prevent these cells from changing again this iteration
        m_lamda[v_ind,] <- 0
      } # n_fromItoJ > 0
    }   # j_u
  }     # i_u
  return(v_U_tm1)
}

# define function to test that spatial output in dt 
# matches input in matrix m_B from a_B posterior
test_matrix_matches <- function(m_B, dt,
  v_times, # = 2016:2020, # unique(blag_pred$df_D$time) 
  v_times_totest # = c(2017, 2018)
  ){

  v_t <- match(v_times_totest, v_times)
  dt[, u_t1 := factor(substr(u_ch, v_t[1], v_t[1]), levels = 1:n_u)]
  dt[, u_t2 := factor(substr(u_ch, v_t[2], v_t[2]), levels = 1:n_u)]
  dt[, u_t12 :=  as.factor(paste0(u_t1, u_t2))]
  # table(dt$u_t12)
  dt[u_t1 == u_t2, u_t12 := NA]
  dt[is.na(u_ch), u_t12 := NA]
   levels(dt$u_t12)
  dt$u_t12 <- droplevels(dt$u_t12)
  v_dt <- table(dt$u_t12)

  v_m <- as.vector(t(m_B))
  v_m <- v_m[v_m > 0]
  #stopifnot(sum(m_B) == sum(table(dt$u_t12)))
  
  return(identical(v_m, as.numeric(v_dt)))
}


# 1. Read in posterior Beta matrix for appropriate region
# n_iter /nchains - burnin / thin posterior samples
# 16749
# (600000-12000)/ 3  / 10 ?
if (region == "en") blag_pred <- tar_read(c_post_B_en)
if (region == "sc") blag_pred <- tar_read(c_post_B_sc)
if (region == "wa") blag_pred <- tar_read(c_post_B_wa)
if (region == "ni") blag_pred <- tar_read(c_post_B_ni)
if (region == "uk") blag_pred <- tar_read(c_post_B_uk)

print(ls.str(blag_pred), max.level = 0)
v_times <- unique(blag_pred$df_D$time) # 1990:2019
n_t <- length(v_times)

a_B_post <- blag_pred$a_B_post
dim(a_B_post)
str(a_B_post)
# set first sample to be MAP
a_B_post[,,, 1] <- blag_pred$a_B_post_map
# set 2nd, 3rd, & 4th sample to be quantiles 0.025, 0.5 , 0.975
a_B_post[,,, 2] <- blag_pred$a_B_post_q[1,,]
a_B_post[,,, 3] <- blag_pred$a_B_post_q[2,,]
a_B_post[,,, 4] <- blag_pred$a_B_post_q[3,,]

# third dimension is time, and has to match likelihood time dimension
# convert area in km^2 to number of grid cells
# units must be km^2
# make function of numeric res
# how many grid cells = 1 km2 at different resolutions
n_cells_per_km2 <- 1e6 / res^2
a_B_post <- round(a_B_post * n_cells_per_km2, 0)
# a_B_post[,, 71, 1:4]
# a_B_post[,, 2, 1]

# 2. check annual likelihood files exist
v_fnames <- here(paste0("data/LikU/", res, "m"), paste0("dt_L_", v_times, ".qs"))
file.exists(v_fnames)

# 3. read current land use and age
fname <- here("data/life_tables", paste0("dt_u_age_", res, "m.qs"))
dt <- qread(fname)
dim(dt)
#str(dt)

# 4. mask dt to land and add region variable
fname <- here("data-raw/mask", paste0("dt_land_", res, "m.qs"))
dt_mask <- qread(fname)
dim(dt_mask)
dim(dt)
#summary(dt)
dt <- data.table(region = dt_mask$country, dt)
#dt[is.na(region), c("u", "age_rev") := NA]
dt[region %!in% i_region, c("u", "age_rev") := NA]
unique(dt$region)
unique(dt_mask$country)
dt[!is.na(u)]

dt[, u_prev := u] # initialise to be the same
dt[, u_ch   := u] # initialise to be the same
#table(dt$age_rev)
#summary(dt)

# 5. Read likelihood from inverse distance to u (Lid)
fname <- here("data/life_tables", paste0("dt_Lid_", res, "m.qs"))
dt_Lid <- qread(fname)
dt_Lid <- dt_Lid * res / 10
# dim(dt_Lid)
# options(digits=10)
# apply(dt_Lid, 2, min, na.rm = TRUE)
# apply(dt_Lid, 2, summary, digits = 10)
# apply(dt_L, 2, summary)
# dt_Lid[!is.na(Lid_1)]

# r <- getRasterTemplate(res = res, crs = crs_OSGB)
# r <- setValues(r, dt_L[, woods])
# plot(r)

# 6. read life tables
a_lamda <- readRDS(here("data/life_tables", "IACS_LifeTables.rds"))


start_time <- Sys.time()
for (i_t in n_t:(2)){
#i_t = n_t - 2
  # 1. read L_static likelihoods for previous year
  dt_L <- qread(v_fnames[(i_t - 1)])

  # 2. lookup life table
  dt[, age_capped := age_rev]
  dt[!(age_rev >= 1 & age_rev <= 10), age_capped := 10]
  dt[, L_1 := a_lamda[cbind(age_capped, u, 1)]]
  dt[, L_2 := a_lamda[cbind(age_capped, u, 2)]]
  dt[, L_3 := a_lamda[cbind(age_capped, u, 3)]]
  dt[, L_4 := a_lamda[cbind(age_capped, u, 4)]]
  dt[, L_5 := a_lamda[cbind(age_capped, u, 5)]]
  dt[, L_6 := a_lamda[cbind(age_capped, u, 6)]]

  # 3. update L1-6 = L_stat * L_dyn
  dt_L <- dt[, .(L_1, L_2, L_3, L_4, L_5, L_6)] *  dt_L
  
  # 4. add inverse-distance likelihood, so min is always > 0
  # apply(dt_L, 2, summary)
  dt_L <- dt_L + dt_Lid
  # apply(dt_L, 2, summary)

  # set all NAs to zero
  dt_L[, 1:n_u][is.na(dt_L[, 1:n_u])] <- 0
  
  # 4. get Beta matrix for this year and use only i_sample (MAP is 1st one)
  m_B <- a_B_post[,, i_t, i_sample]
  diag(m_B) <- 0
  # sum(m_B)

  # 5. sample the previous year's land use
  dt[, u_prev := getUt(u, m_B, dt_L)]
  #dt[!is.na(u)]
  
  # 6. prepend it to the character string
  dt[, u_ch := stri_c(u_prev, u_ch, sep = "")]

  # 7. update age_rev
  dt[u != u_prev, age_rev := 1]
  dt[u == u_prev, age_rev := age_rev + 1]
  
  # 8. and finally reset u to previous land use (because we're going backwards)
  dt[, u := u_prev]  ##* WIP TEMP comment out
   
  # some checks
  # dim(dt[u_prev != u])
  # dt[u_prev == u]
  # ind <- which(v_U_tm1 != v_U_t)
  # dt[ind, u]
  # dt[ind, u_tm1]
  # v_U_t[ind]
  # v_U_tm1[ind]
  # dt$v_U_tm1[ind]
} # n_t
print(paste0("time elapsed: ", Sys.time() - start_time))

dt[!is.na(u)]

# verify output
v_times_post <- 1950:2020 # in a_B_post array
for (time_totest in 1951:2020){ # start at dt start+1
#time_totest = 2020
  v_t <- match(time_totest, v_times_post)
  m_B <- a_B_post[,, v_t, i_sample]
  diag(m_B) <- 0
  print(test_matrix_matches(m_B, dt, 
    v_times = 1950:2020, # in dt u_ch
    v_times_totest = c(time_totest - 1, time_totest)))
}

# save output
fname <- here(dir_output, paste0("dt_u_", region, "_smp", i_sample, "_", res, "m.qs"))
qsave(dt[, .(u_ch)], file = fname)

# # re-read file in if needed
# fname <- here(dir_output, paste0("dt_u_", region, "_smp", i_sample, "_", res, "m.qs"))
# dt <- qread(file = fname)

dt[!is.na(u_ch)]
dt[!is.na(u_t12)]
nchar(dt[!is.na(u_ch), u_ch][1])

# condense into unique land-use vector format
dt_luv <- dt[, .N, by = u_ch]
dt_luv[, area := N * res^2 / 1e6] # res in m^2, area in km^2
dt_luv[, start_time := v_times[1]]
dt_luv[, end_time := v_times[length(v_times)]]
dt_luv <- dt_luv[order(-area)]
dt_luv <- dt_luv[!is.na(u_ch)] # remove NA, only first row with all non-land grid cells
sum(dt_luv[, area])
dim(dt_luv)
# save output
fname <- here(dir_output, paste0("dt_luv_", region, "_smp", i_sample, "_", res, "m.qs"))
qsave(dt_luv, file = fname)

# dt_luv <- qread(file = fname)
library(stringr)
# remove the six vectors with no change
ind_unchanging <- c(
  str_which(dt_luv$u_ch, "1{71}"),
  str_which(dt_luv$u_ch, "2{71}"),
  str_which(dt_luv$u_ch, "3{71}"),
  str_which(dt_luv$u_ch, "4{71}"),
  str_which(dt_luv$u_ch, "5{71}"),
  str_which(dt_luv$u_ch, "6{71}")
)
dt_luv <- dt_luv[-ind_unchanging]

# identify to crop-grass rotational vectors
ind_cgc <- str_which(dt_luv$u_ch, "2+3+2+")
ind_gcg <- str_which(dt_luv$u_ch, "3+2+3+")

dt_luv$id <- as.integer(rownames(dt_luv))

dt_luv_wide <- dt_luv %>% separate(u_ch, into = as.character(v_times), sep = 1:(n_t-1), convert = TRUE, remove = FALSE)
names(dt_luv_wide)
dt_luv_long <- 
  pivot_longer(dt_luv_wide, cols = as.character(v_times), names_to = "time", values_to = "u")
names(dt_luv_long)
dt_luv_long$time <- as.integer(dt_luv_long$time)


p <- ggplot(subset(dt_luv_long, id < 100), aes(time, u))
#p <- ggplot(subset(dt_luv_long, id %in% ind_cgc & area > 5), aes(time, u))
p <- p + geom_path(aes(size = area, alpha = area, colour = as.factor(id)))
p <- p + scale_y_continuous(breaks = 1:6, labels = names_u)
p <- p + xlab("Year") + ylab("Land use")
p <- p + guides(colour=FALSE, alpha=FALSE, size=FALSE)
p_v1 <- p

p <- ggplot(subset(dt_luv_long, id >= 1 & id <= 42), aes(time, u))
#p <- ggplot(subset(dt_luv_long, id %in% ind_cgc & area > 7), aes(time, u))
#p <- p + geom_path(aes(size = area_lu_vector, colour = as.factor(id)))
p <- p + geom_path(aes(size = area))
p <- p + scale_y_continuous(breaks = 1:6, labels = names_u)
p <- p + facet_wrap( ~ id, scales = "fixed")
p <- p + guides(colour=FALSE)
p <- p + xlab("Year") + ylab("Land use")
p <- p + scale_x_continuous(breaks=c(1975,2000))
p_v2 <- p

r <- getRasterTemplate(domain = "UK", res = res, crs = crs_OSGB)
r_u_t1   <- setValues(r, dt[, u_t1])
r_u_t2   <- setValues(r, dt[, u_t2])
r_u_t12   <- setValues(r, dt[, u_t12])

dt_plot <- as.data.table(as.data.frame(r_u_t12, xy = TRUE))
dt_plot <- data.table(dt_plot, dt[, .(u_t1, u_t2)])
names(dt_plot)[3] <- "u_u"
dt_plot <- dt_plot[!is.na(u_u)]

#p <- ggplot(data = dt_plot[sample(1:dim(dt_plot)[1], 2000)], aes(x,y)) + 
p_u_d <- ggplot(data = dt_plot, aes(x,y)) + 
  geom_point(aes(colour = u_t1), shape = "square", size = 3) +
  geom_point(aes(colour = u_t2), shape = "circle", size = 1.5) +
  scale_colour_manual(name = "Square = previous \n Circle = new", values = colour_u, labels = names_u) #+
  #coord_equal() +  theme_map() +  theme(legend.position = "right")

v_times_toplot <- c(2019, 2020)
s <- stack(r_u_t1, r_u_t2)
st_U <- st_as_stars(s)
st_U <- st_set_dimensions(st_U, "band", values = v_times_toplot, names = "time")
names(st_U) <- "U"
st_U$U <- factor(st_U$U)
levels(st_U) <- names_u
#plot(st_U, interpolate = TRUE)

p_u_t <- ggplot() + 
  geom_stars(data = st_U[,,, c(1, 2)], 
  downsample = c(1, 1, 1)) + 
  facet_wrap("time") +
  scale_fill_manual(values = colour_u, labels = names_u) +
  coord_equal() +  theme_map() +
  theme(legend.position = "right")

fname <- here(dir_output, paste0("r_uu_", region, "_smp", i_sample, "_", res, "m.pdf"))
pdf(fname)
  p_u_t
  p_u_d
  p_v1
  p_v2
dev.off()

# older graphics
# plot some years
# v_times <- 2016:2020 # unique(blag_pred$df_D$time) 
# v_times_toplot <- c(2017, 2018)
# v_t <- match(v_times_toplot, v_times)
# n_t <- length(v_times)

# st_U <- st_as_stars(r_u_t12)
# st_U <- st_set_dimensions(st_U, "band", values = v_times_toplot, names = "time")
# names(st_U) <- "U"
# st_U$U <- factor(st_U$U)
# levels(st_U) <- names_u
# plot(st_U, interpolate = TRUE)

# p <- ggplot() + 
  # geom_stars(data = st_U[,,, c(1, 2, 3)], 
  # downsample = c(1, 1, 1)) + 
  # facet_wrap("time") +
  # scale_fill_manual(values = colour_u, labels = names_u) +
  # coord_equal() +  theme_map() +
  # theme(legend.position = "right")
# p


#** WIP
# Beta matrix has some odd values rough to urban
# use newer B_post with half-normal prior
# may be oddity of 10 km res?
# check getUt function is sensible - with chess?
# run getBLAG_fromU function on it to check 
