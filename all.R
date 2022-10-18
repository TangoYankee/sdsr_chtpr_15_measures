library(sf)
library(tmap)
library(spdep) |> suppressPackageStartupMessages()
library(parallel)
library(spData)
# To use spDataLarge
## install.packages("raster")
## install.packages("spDataLarge", repos = "https://nowosad.github.io/drat/", type = "source")
library(spDataLarge)
library(rgeoda)

## Chapter 14 setup
data(pol_pres15, package = "spDataLarge")

if (!all(st_is_valid(pol_pres15)))
  pol_pres15 <- st_make_valid(pol_pres15)

pol_pres15 |> poly2nb(queen = TRUE) -> nb_q

(nb_q |> 
    nb2listw(style = "B") -> lw_q_B) |> 
  spweights.constants() |> 
  data.frame() |> 
  subset(select = c(n, S0, S1, S2))

(nb_q |> 
    nb2listw(style = "W") -> lw_q_W) |> 
  spweights.constants() |> 
  data.frame() |> 
  subset(select = c(n, S0, S1, S2))

pol_pres15 |> 
  st_geometry() |> 
  st_centroid(of_largest_polygon = TRUE) -> coords 

(coords |> dnearneigh(0, 18300) -> nb_d183)

nb_d183 |> 
  nbdists(coords) |> 
  lapply(function(x) 1/(x/1000)) -> gwts
(nb_d183 |> nb2listw(glist=gwts, style="B") -> lw_d183_idw_B) |> 
  spweights.constants() |> 
  data.frame() |> 
  subset(select=c(n, S0, S1, S2))

## Start Chapter 15
### 15.1
glance_htest <- function(ht) c(ht$estimate, 
                               "Std deviate" = unname(ht$statistic), 
                               "p.value" = unname(ht$p.value))
set.seed(1)
(pol_pres15 |> 
    nrow() |>
    rnorm() -> x) |> 
  moran.test(lw_q_B, randomisation = FALSE,
             alternative = "two.sided") |> 
  glance_htest()

beta <- 0.0015
coords |> 
  st_coordinates() |> 
  subset(select = 1, drop = TRUE) |> 
  (function(x) x/1000)() -> t
(x + beta * t -> x_t) |> 
  moran.test(lw_q_B, randomisation = FALSE,
             alternative = "two.sided") |> 
  glance_htest()

lm(x_t ~ t) |> 
  lm.morantest(lw_q_B, alternative = "two.sided") |> 
  glance_htest()

# 15.2
(pol_pres15 |> 
    st_drop_geometry() |> 
    subset(select = types, drop = TRUE) -> Types) |> 
  table()

Types |> joincount.multi(listw = lw_q_B)

Types |> joincount.multi(listw = lw_d183_idw_B)

pol_pres15 |> 
  st_drop_geometry() |> 
  subset(select = I_turnout, drop = TRUE) -> I_turnout

I_turnout |> moran.test(listw = lw_q_B, randomisation = FALSE) |> 
  glance_htest()

lm(I_turnout ~ 1, pol_pres15) |> 
  lm.morantest(listw = lw_q_B) |> 
  glance_htest()

(I_turnout |> 
    moran.test(listw = lw_q_B) -> mtr) |> 
  glance_htest()

set.seed(1)
I_turnout |> 
  moran.mc(listw = lw_q_B, nsim = 999, 
           return_boot = TRUE) -> mmc

c("Permutation bootstrap" = var(mmc$t), 
  "Analytical randomisation" = unname(mtr$estimate[3]))

## 15.3
I_turnout |> 
  moran.plot(listw = lw_q_W, labels = pol_pres15$TERYT, 
             cex = 1, pch = ".", xlab = "I round turnout", 
             ylab = "lagged turnout") -> infl_W

pol_pres15$hat_value <- infl_W$hat
tm_shape(pol_pres15) + tm_fill("hat_value")

I_turnout |> 
  localmoran(listw = lw_q_W) -> locm

all.equal(sum(locm[,1])/Szero(lw_q_W), 
          unname(moran.test(I_turnout, lw_q_W)$estimate[1]))

pva <- function(pv) cbind("none" = pv, 
                          "FDR" = p.adjust(pv, "fdr"), "BY" = p.adjust(pv, "BY"),
                          "Bonferroni" = p.adjust(pv, "bonferroni"))
locm |> 
  subset(select = "Pr(z != E(Ii))", drop = TRUE) |> 
  pva() -> pvsp
f <- function(x) sum(x < 0.005)
apply(pvsp, 2, f)

invisible(spdep::set.coresOption(max(detectCores()-1L, 1L)))
I_turnout |> 
  localmoran_perm(listw = lw_q_W, nsim = 9999, 
                  iseed = 1) -> locm_p

locm_p |> 
  subset(select = "Pr(z != E(Ii))", drop = TRUE) |> 
  pva() -> pvsp
apply(pvsp, 2, f)

locm_p |> 
  subset(select = "Pr(z != E(Ii)) Sim", drop = TRUE) |> 
  pva() -> pvsp
apply(pvsp, 2, f)

pol_pres15$locm_pv <- p.adjust(locm[, "Pr(z != E(Ii))"], "fdr")
pol_pres15$locm_std_pv <- p.adjust(locm_p[, "Pr(z != E(Ii))"], 
                                   "fdr")
pol_pres15$locm_p_pv <- p.adjust(locm_p[, "Pr(z != E(Ii)) Sim"],
                                 "fdr")

tm_shape(pol_pres15) +
  tm_fill(c("locm_pv", "locm_std_pv", "locm_p_pv"), 
          breaks=c(0, 0.0005, 0.001, 0.005, 0.01, 
                   0.05, 0.1, 0.2, 0.5, 0.75, 1), 
          title = "Pseudo p-values\nLocal Moran's I",
          palette="-YlOrBr") +
  tm_facets(free.scales = FALSE, ncol = 2) +
  tm_layout(panel.labels = c("Analytical conditional",
                             "Permutation std. dev.",
                             "Permutation rank"))

quadr <- attr(locm, "quadr")$mean
a <- table(addNA(quadr))
locm |> hotspot(Prname="Pr(z != E(Ii))", cutoff = 0.005, 
                droplevels=FALSE) -> pol_pres15$hs_an_q
locm_p |> hotspot(Prname="Pr(z != E(Ii))", cutoff = 0.005, 
                  droplevels=FALSE) -> pol_pres15$hs_ac_q 
locm_p |> hotspot(Prname="Pr(z != E(Ii)) Sim", cutoff = 0.005,
                  droplevels = FALSE) -> pol_pres15$hs_cp_q
b <- table(addNA(pol_pres15$hs_an_q))
c <- table(addNA(pol_pres15$hs_ac_q))
d <- table(addNA(pol_pres15$hs_cp_q))
t(rbind("Moran plot quadrants" = a, "Analytical cond." = b, 
        "Permutation std. cond." = c, "Permutation rank cond." = d))
#           Moran plot quadrants Analytical cond.
# Low-Low                   1040               53
# High-Low                   264                0
# Low-High                   213                0
# High-High                  978               96
# <NA>                         0             2346
#           Permutation std. cond. Permutation rank cond.
# Low-Low                       53                     55
# High-Low                       0                      0
# Low-High                       0                      0
# High-High                     96                     70
# <NA>                        2346                   2370

pol_pres15$hs_an_q <- droplevels(pol_pres15$hs_an_q)
pol_pres15$hs_ac_q <- droplevels(pol_pres15$hs_ac_q)
pol_pres15$hs_cp_q <- droplevels(pol_pres15$hs_cp_q)

### Figure 15.4
tm_shape(pol_pres15) +
  tm_fill(c("hs_an_q", "hs_ac_q", "hs_cp_q"),
          colorNA = "grey95", textNA="Not \"interesting\"",
          title = "Turnout hotspot status\nLocal Moran's I",
          palette = RColorBrewer::brewer.pal(4, "Set3")[-c(2,3)]) +
  tm_facets(free.scales = FALSE, ncol = 2) +
  tm_layout(panel.labels = c("Analytical conditional",
                             "Permutation std. cond.",
                             "Permutation rank cond."))

lm(I_turnout ~ 1) -> lm_null
lm_null |> localmoran.sad(nb = nb_q, style = "W",
                          alternative = "two.sided") |>
  summary() -> locm_sad_null

lm(I_turnout ~ 1, weights = pol_pres15$I_entitled_to_vote) ->
  lm_null_weights
lm_null_weights |>
  localmoran.sad(nb = nb_q, style = "W",
                 alternative = "two.sided") |>
  summary() -> locm_sad_null_weights

lm(I_turnout ~ Types, weights=pol_pres15$I_entitled_to_vote) ->
  lm_types
lm_types |> localmoran.sad(nb = nb_q, style = "W",
                           alternative = "two.sided") |>
  summary() -> locm_sad_types

locm_sad_null |> hotspot(Prname="Pr. (Sad)",
                         cutoff=0.005) -> pol_pres15$locm_sad0
locm_sad_null_weights |> hotspot(Prname="Pr. (Sad)",
                                 cutoff = 0.005) -> pol_pres15$locm_sad1
locm_sad_types |> hotspot(Prname="Pr. (Sad)",
                          cutoff = 0.005) -> pol_pres15$locm_sad2

### Figure 15.5
tm_shape(pol_pres15) +
  tm_fill(c("hs_cp_q", "locm_sad0", "locm_sad1",  "locm_sad2"),
          colorNA = "grey95", textNA = "Not \"interesting\"",
          title = "Turnout hotspot status\nLocal Moran's I",
          palette =RColorBrewer::brewer.pal(4, "Set3")[c(1, 4, 2)]) +
  tm_facets(free.scales = FALSE, ncol = 2) + 
  tm_layout(panel.labels = c("Permutation rank", 
                             "saddlepoint null", "saddlepoint weighted null", 
                             "saddlepoint weighted types"))

rbind(null = append(table(addNA(pol_pres15$locm_sad0)),
                    c("Low-High" = 0), 1),
      weighted = append(table(addNA(pol_pres15$locm_sad1)),
                        c("Low-High" = 0), 1),
      type_weighted = append(table(addNA(pol_pres15$locm_sad2)),
                             c("Low-High" = 0), 1))

lm_types |> localmoran.exact(nb = nb_q, style = "W", 
                             alternative = "two.sided", useTP=TRUE, truncErr=1e-8) |> 
  as.data.frame() -> locm_ex_types

locm_ex_types |> hotspot(Prname = "Pr. (exact)",
                         cutoff = 0.005) -> pol_pres15$locm_ex

### Figure 15.6
tm_shape(pol_pres15) +
  tm_fill(c("locm_sad2", "locm_ex"), colorNA = "grey95",
          textNA = "Not \"interesting\"", 
          title = "Turnout hotspot status\nLocal Moran's I",
          palette = RColorBrewer::brewer.pal(4, "Set3")[c(1, 4, 2)]) +
  tm_facets(free.scales = FALSE, ncol = 2) +
  tm_layout(panel.labels = c("saddlepoint weighted types",
                             "Exact weighted types"))

table(Saddlepoint = addNA(pol_pres15$locm_sad2),
      exact = addNA(pol_pres15$locm_ex))

I_turnout |> 
  localG(lw_q_W, return_internals = TRUE) -> locG

I_turnout |> 
  localG_perm(lw_q_W, nsim = 9999, iseed = 1) -> locG_p

cor(cbind(localG=attr(locG, "internals")[, "Pr(z != E(Gi))"], 
          attr(locG_p, "internals")[, c("Pr(z != E(Gi))", 
                                        "Pr(z != E(Gi)) Sim")]))

I_turnout |> 
  localC_perm(lw_q_W, nsim=9999, iseed=1) -> locC_p

cor(attr(locC_p, "pseudo-p")[, c("Pr(z != E(Ci))",
                                 "Pr(z != E(Ci)) Sim")])

locC_p |> hotspot(Prname = "Pr(z != E(Ci)) Sim",
                  cutoff = 0.005) -> pol_pres15$hs_C
locG_p |> hotspot(Prname = "Pr(z != E(Gi)) Sim",
                  cutoff = 0.005) -> pol_pres15$hs_G

#### Figure 15.7
m1 <- tm_shape(pol_pres15) +
  tm_fill("hs_cp_q", 
          palette = RColorBrewer::brewer.pal(4, "Set3")[-c(2,3)],
          colorNA = "grey95", textNA = "Not \"interesting\"",
          title = "Turnout hotspot status\nLocal Moran I") + 
  tm_layout(legend.outside=TRUE, legend.outside.position="bottom")
m2 <- tm_shape(pol_pres15) +
  tm_fill("hs_G",
          palette = RColorBrewer::brewer.pal(4, "Set3")[-c(2,3)],
          colorNA = "grey95", textNA="Not \"interesting\"",
          title = "Turnout hotspot status\nLocal Getis/Ord G") +
  tm_layout(legend.outside=TRUE, legend.outside.position="bottom")
m3 <- tm_shape(pol_pres15) +
  tm_fill("hs_C",
          palette = RColorBrewer::brewer.pal(4, "Set3")[c(4, 1, 3)],
          colorNA = "grey95", textNA="Not \"interesting\"",
          title = "Turnout hotspot status\nLocal Geary C") +
  tm_layout(legend.outside=TRUE, legend.outside.position="bottom")
tmap_arrange(m1, m2, m3, nrow=1)

pol_pres15 |> 
  st_drop_geometry() |> 
  subset(select = II_turnout) |> 
  localC_perm(lw_q_W, nsim=9999, iseed=1) -> locC_p_II

locC_p_II |> hotspot(Prname = "Pr(z != E(Ci)) Sim",
                     cutoff = 0.005) -> pol_pres15$hs_C_II

pol_pres15 |> 
  st_drop_geometry() |> 
  subset(select = c(I_turnout, II_turnout)) |>
  localC_perm(lw_q_W, nsim=9999, iseed=1) -> locMvC_p

all.equal(locMvC_p, (locC_p+locC_p_II)/2,
          check.attributes = FALSE)

locMvC_p |> hotspot(Prname = "Pr(z != E(Ci)) Sim",
                    cutoff = 0.005) -> pol_pres15$hs_MvC

#### Figure 15.8
m3 <- tm_shape(pol_pres15) +
  tm_fill("hs_C", 
          palette = RColorBrewer::brewer.pal(4, "Set3")[c(4, 1, 3, 2)],
          colorNA = "grey95", textNA = "Not \"interesting\"",
          title = "First round turnout\nLocal Geary C") +
  tm_layout(legend.outside=TRUE, legend.outside.position="bottom")
m4 <- tm_shape(pol_pres15) +
  tm_fill("hs_C_II", 
          palette = RColorBrewer::brewer.pal(4, "Set3")[c(4, 1, 3, 2)], 
          colorNA = "grey95", textNA = "Not \"interesting\"",
          title="Second round turnout\nLocal Geary C") +
  tm_layout(legend.outside=TRUE, legend.outside.position="bottom")
m5 <- tm_shape(pol_pres15) +
  tm_fill("hs_MvC", 
          palette = RColorBrewer::brewer.pal(4, "Set3")[c(4, 1)],
          colorNA = "grey95", textNA = "Not \"interesting\"",
          title = "Both rounds turnout\nLocal Multivariate Geary C") +
  tm_layout(legend.outside=TRUE, legend.outside.position="bottom")
tmap_arrange(m3, m4, m5, nrow=1)

table(droplevels(interaction(addNA(pol_pres15$hs_C),
                             addNA(pol_pres15$hs_C_II), sep=":")), 
      addNA(pol_pres15$hs_MvC))

## The rgeoda package
Geoda_w <- queen_weights(pol_pres15)
summary(Geoda_w)

lisa <- local_multigeary(Geoda_w, 
                         pol_pres15[c("I_turnout", "II_turnout")], 
                         cpu_threads = max(detectCores() - 1, 1),
                         permutations = 99999, seed = 1)

all.equal(card(nb_q), lisa_num_nbrs(lisa), 
          check.attributes = FALSE)
all.equal(lisa_values(lisa), c(locMvC_p),
          check.attributes = FALSE)

apply(attr(locMvC_p, "pseudo-p")[,c("Pr(z != E(Ci)) Sim", 
                                    "Pr(folded) Sim")], 2, range)

locMvC_p |> hotspot(Prname = "Pr(folded) Sim",
                    cutoff = 0.0025) -> pol_pres15$hs_MvCa

mvc <- factor(lisa_clusters(lisa), levels=0:2,
              labels = lisa_labels(lisa)[1:3])
is.na(mvc) <- p.adjust(lisa_pvalues(lisa), "fdr") >= 0.0025
pol_pres15$geoda_mvc <- droplevels(mvc)

addmargins(table(spdep = addNA(pol_pres15$hs_MvCa),
                 rgeoda = addNA(pol_pres15$geoda_mvc)))

#### 15.9
m5 <- tm_shape(pol_pres15) +
  tm_fill("hs_MvCa", 
          palette = RColorBrewer::brewer.pal(4, "Set3")[c(4, 1)],
          colorNA = "grey95", textNA = "Not \"interesting\"",
          title = "Both rounds turnout (spdep)\nLocal Multivariate Geary C")
m6 <- tm_shape(pol_pres15) +
  tm_fill("geoda_mvc", 
          palette = RColorBrewer::brewer.pal(4, "Set3")[c(4, 1)],
          colorNA = "grey95", textNA = "Not \"interesting\"",
          title="Both rounds turnout (rgeoda)\nLocal Multivariate Geary C")
tmap_arrange(m5, m6, nrow=1)
