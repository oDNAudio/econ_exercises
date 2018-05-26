library(BMS)
library(ggplot2)
library(ggthemes)

source("exercise_3_ssvs.R")

data(datafls)

Y = datafls[, 1]
X = datafls[, 2:ncol(datafls)]
rm(datafls)

tau0s = c(1, 1e-2, 1e-5, 1e-15)
tau1_scales = c(100, 1000)

pips = vector("list", length(tau0s) * length(tau1_scales))

i = 1
for(scale in tau1_scales) {
  for(tau in tau0s) {
    pips[[i]] = ssvs(Y, X, tau0 = tau, tau1 = tau * scale, standardise = FALSE)[[1]]
    names(pips)[i] = paste0(tau, " & *", scale)
    i = i + 1
  }
}

pips_z = vector("list", length(tau0s) * length(tau1_scales))

i = 1
for(scale in tau1_scales) {
  for(tau in tau0s) {
    pips_z[[i]] = ssvs(Y, X, tau0 = tau, tau1 = tau * scale, standardise = TRUE)[[1]]
    names(pips_z)[i] = paste0(tau, " & *", scale)
    i = i + 1
  }
}

pips_auto = vector("list", 1)
pips_auto = ssvs(Y, X)[[1]]

# Demonstrate graphically

df = data.frame(pips_z[[2]], pips_z[[4]], pips_z[[6]])
df2 = data.frame(pips[[2]], pips[[4]], pips[[6]])
names(df) = names(df2) = c("0.01(*100)", "1e-15(*100)", "0.01(*1000)")
df$id = df2$id = as.numeric(rownames(df))

p1 = ggplot(df, aes(x = id, y = `0.01(*100)`)) +
  geom_point(aes(colour = "z, 0.01(*100)")) +
  geom_smooth(aes(colour = "z, 0.01(*100)"), alpha = 0.2) +
  geom_point(aes(y = `1e-15(*100)`, colour = "z, 1e-15(*100)")) +
  geom_smooth(aes(y = `1e-15(*100)`, colour = "z, 1e-15(*100)"), alpha = 0.2) +
  geom_point(aes(y = `0.01(*1000)`, colour = "z, 0.01(*1000)")) +
  geom_smooth(aes(y = `0.01(*1000)`, colour = "z, 0.01(*1000)"), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 1), expand = 0) +
  theme_fivethirtyeight() +
  scale_color_gdocs(name = "")
p1

p2 = ggplot(df2, aes(x = id, y = `0.01(*100)`)) +
  geom_point(aes(colour = "no-z, 0.01(*100)")) +
  geom_smooth(aes(colour = "no-z, 0.01(*100)"), alpha = 0.2) +
  geom_point(aes(y = `1e-15(*100)`, colour = "no-z, 1e-15(*100)")) +
  geom_smooth(aes(y = `1e-15(*100)`, colour = "no-z, 1e-15(*100)"), alpha = 0.2) +
  geom_point(aes(y = `0.01(*1000)`, colour = "no-z, 0.01(*1000)")) +
  geom_smooth(aes(y = `0.01(*1000)`, colour = "no-z, 0.01(*1000)"), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 1), expand = 0) +
  theme_fivethirtyeight() +
  scale_color_gdocs(name = "")
p2

df3 = data.frame(pips_auto, pips_z[[2]], pips_z[[6]])
names(df3) = c("automatic", "0.01(*100)", "0.01(*1000)")
df3$id = as.numeric(rownames(df3))

p3 = ggplot(df3, aes(x = id, y = automatic)) +
  geom_point(aes(colour = "automatic")) +
  geom_smooth(aes(colour = "automatic"), alpha = 0.2) +
  geom_point(aes(y = `0.01(*100)`, colour = "z, 0.01(*100)")) +
  geom_smooth(aes(y = `0.01(*100)`, colour = "z, 0.01(*100)"), alpha = 0.2) +
  geom_point(aes(y = `0.01(*1000)`, colour = "z, 0.01(*1000)")) +
  geom_smooth(aes(y = `0.01(*1000)`, colour = "z, 0.01(*1000)"), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 1), expand = 0) +
  theme_fivethirtyeight() +
  scale_color_gdocs(name = "")
p3

df4 = data.frame(pips_auto, pips_z[[2]], pips[[2]])
names(df4) = c("automatic", "standardised", "neither")
df4$id = as.numeric(rownames(df4))

p4 = ggplot(df4, aes(x = id, y = automatic)) +
  geom_point(aes(colour = "automatic")) +
  geom_smooth(aes(colour = "automatic"), alpha = 0.2) +
  geom_point(aes(y = standardised, colour = "z, 0.01(*100)")) +
  geom_smooth(aes(y = standardised, colour = "z, 0.01(*100)"), alpha = 0.2) +
  geom_point(aes(y = neither, colour = "no-z, 0.01(*100)")) +
  geom_smooth(aes(y = neither, colour = "no-z, 0.01(*100)"), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 1), expand = 0) +
  theme_fivethirtyeight() +
  scale_color_gdocs(name = "")
p4

plot_grid(p1, p2, p3, p4)
