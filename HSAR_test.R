pacman::p_load(tidyverse, spdep, spatialreg, tmap, RColorBrewer,
               mapview, viridisLite, leaflet.extras2)

rm(list = ls())

load("twenty_two.RData")

# url <- "https://cran.r-project.org/src/contrib/Archive/HSAR/HSAR_0.5.tar.gz"
#install.packages(url, type="source", repos=NULL)
library(HSAR)

# Define the random effect matrix
model.data <- twenty_two[order(twenty_two$MSOA11CD),]
head(model.data,50)

# the number of LSOAS within each MSOAS
MM <- as.data.frame(table(model.data$MSOA11CD))
# the total number of neighbourhood, 100
Utotal <- dim(MM)[1]
Unum <- MM[,2]
Uid <- rep(c(1:Utotal),Unum)

n <- nrow(model.data)
Delta <- matrix(0,nrow=n,ncol=Utotal)
for(i in 1:Utotal) {
  Delta[Uid==i,i] <- 1
}
rm(i)
# Delta[1:50,1:10]
Delta <- as(Delta,"dgCMatrix")

msoas <- read_sf("MSOA_2011_EW_BGC_V3.shp")
sheff_msoas <- twenty_two$MSOA11CD
msoas <- msoas %>% 
  filter(MSOA11CD %in% sheff_msoas)

# extract the MSOA level spatial weights matrix using the queen's rule
nb.list <- spdep::poly2nb(msoas)
mat.list <- spdep::nb2mat(nb.list,style="W")
M <- as(mat.list,"dgCMatrix")

# extract the LSOA level spatial weights matrix
lsoa_cent <- st_centroid(twenty_two)
nb.25 <- spdep::dnearneigh(lsoa_cent,0,6000)
# to a weights matrix
dist.25 <- spdep::nbdists(nb.25,lsoa_cent)
dist.25 <- lapply(dist.25,function(x) exp(-0.5 * (x / 2500)^2))
mat.25 <- spdep::nb2mat(nb.25,glist=dist.25,style="W")
W <- as(mat.25,"dgCMatrix")

# imd scores
imd <- read_csv("imd_scores.csv")
imd <- imd %>% select(lsoa_code, imd_score)

model.data <- model.data %>% 
  left_join(imd, by = c("LSOA11CD" = "lsoa_code"))
twenty_two <- twenty_two %>% 
  left_join(imd, by = c("LSOA11CD" = "lsoa_code"))

## run the hsar() function
res.formula <- price ~ imd_score

betas= coef(lm(formula=res.formula,data=twenty_two))
pars=list( rho = 0.5,lambda = 0.5, sigma2e = 2.0, sigma2u = 2.0, betas = betas )

# }
# NOT RUN {
res <- hsar(res.formula,data=model.data,W=W,M=M,Delta=Delta,
            burnin=500, Nsim=1000, thinning = 1, parameters.start=pars)
summary(res)

# visualise the district level random effect
library(classInt)
x <- as.numeric(res$Mus)
break_n <- 7
breaks <- classIntervals(x,break_n,"equal")$brks
groups <- cut(x,breaks,include.lowest=TRUE,labels=FALSE)
palette <- viridisLite::viridis(option = "turbo",n = break_n, alpha = .5)
# plot(msoas["MSOA11CD"],col=palette[groups],border="grey")
# }
# NOT RUN {
# }

tmap_mode("view")
msoas$groups <- cut(x,breaks,include.lowest=TRUE,labels=FALSE)
msoas$x <- x

tm_shape(msoas) + 
  tm_fill(col = "groups", title = "",
          palette =  palette) +
  tm_legend(text.size = 1)  + tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE,  title = "HSAR Clusters")  +
  tm_layout(legend.outside = TRUE)

res$R_Squared
res$DIC

# mlm model ---------------------------------------------------------------

pacman::p_load(lme4, lmerTest, jtools)

mlm_mod <- lmer(price ~ imd_score + (1|MSOA11CD),
                data = model.data, REML = FALSE)

summ(mlm_mod)

ran_ints <- tibble(ranefs = ranef(mlm_mod)$MSOA11CD$`(Intercept)`,
                   MSOA11CD = row.names(ranef(mlm_mod)$MSOA11CD))

msoas <- msoas %>% 
  left_join(ran_ints, by = "MSOA11CD")

breaks2 <- classIntervals(msoas$ranefs,break_n,"equal")$brks
groups2 <- cut(msoas$ranefs,breaks2,include.lowest=TRUE,labels=FALSE)

msoas$groups2 <- cut(msoas$ranefs,breaks2,include.lowest=TRUE,labels=FALSE)

tm_shape(msoas) + 
  tm_fill(col = "groups2", title = "",
          palette =  palette) +
  tm_legend(text.size = 1)  + tm_borders(alpha = 0.5) +
  tm_layout(frame = FALSE,  title = "MLM Clusters")  +
  tm_layout(legend.outside = TRUE)

# ggplot ----------------------------------------------------------------

pacman::p_load(patchwork)

p1 <- msoas %>% 
  ggplot(aes(fill = scale(x))) +
  geom_sf(alpha = 0.8, show.legend = FALSE) +
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  labs(fill = "Random\neffects")

p2 <- msoas %>% 
  ggplot(aes(fill = scale(ranefs))) +
  geom_sf(alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo") +
  theme_void() +
  labs(fill = "Random\neffects")

p1 + p2

# leaflet slider map ------------------------------------------------------

msoas_l <- st_transform(msoas, 4326)

Sheffield <- msoas_l %>% 
  mutate(HSAR = scale(x)[,1],
         MLM = scale(ranefs)[,1])

pal <- colorRampPalette(viridis(n=8))
at1 <- seq(min(Sheffield$HSAR),max(Sheffield$HSAR),length.out = 7)
at2 <- seq(min(Sheffield$MLM),max(Sheffield$MLM),length.out = 7)

m1 <- mapview(Sheffield, zcol = "HSAR", map.types = "CartoDB.Positron",
              col.regions = pal, at = at1)
m2 <- mapview(Sheffield, zcol = "MLM", map.types = "CartoDB.Positron",
              col.regions = pal, at = at2)

m1 | m2

