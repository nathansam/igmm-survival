### R code to create the figures used in 03/08/2020 IGMM statistical seminar ###
### series presentation on survival analysis by Nathan Constantine-Cooke     ###

################################################################################
################################# PREAMBLE #####################################
################################################################################

library(ggplot2)
library(survival) # Main R survival analysis package
library(survminer) # Used for survival analysis visualisations

# Weibull hazard and survival functions
weibull_survival <- function(x, p, lambda) {
        exp(-lambda * (t ^ p))
}

weibull_hazard <- function(t, p, lambda) {
        lambda * p * t ^ (p - 1)
}

# ggplot2 theme for plots
custom_theme <- theme_classic() + theme(axis.text = element_text(size = 14),
                                        axis.title = element_text(size = 18,
                                                                  face = "bold"))
# Time values used for all continuous (non step-wise) plots
t <- seq(0, 5, by = 0.01)

# Create figures directory
if (dir.exists("Figures") == FALSE) {
        dir.create("Figures")
}

################################################################################
#################################### FIGURES ###################################
################################################################################

## Theoretical survival plot ##
st <- weibull_survival(t, 2.4, 0.4)
plot.dat <- data.frame(t, st)
p <- ggplot(aes(x = t, y = st), data = plot.dat)
p <- p + geom_line()
p <- p + custom_theme
p <- p + ylab("S(t)")
ggsave("Figures/Survival.svg",
       p,
       device = "svg",
       width = 15.5 / 1.5,
       height = 10.3 / 1.5)


## "real-world" survival plot ##
p <- ggsurvplot(fit = survfit(Surv(time, status) ~ 1, data = lung),
                xlab = "Days",
                ylab = expression(hat(S)(t)),
                conf.int = FALSE,
                ggtheme = custom_theme,
                palette = "black")
# Print because can't pass p directly to ggsave()
print(p)
ggsave("Figures/Survival2.svg",
       device = "svg",
       width = 15.5 / 1.5,
       height = 10.3 / 1.5)

## Stratified Kaplan-Meier plot ##
data("ovarian")
surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
ovarian$rx <- factor(ovarian$rx, levels = c(1, 2), labels = c("A", "B"))
fit1 <- survfit(surv_object ~ rx, data = ovarian)
ggsurvplot(fit1, data = ovarian, pval = FALSE, conf.int = FALSE)
ggsave("Figures/Kaplan-Meier.svg",
       device = "svg",
       width = 15.5 / 1.5,
       height = 10.3 / 1.5)


## Weibull survival plots w/ varied p values ##
for (i in c(0.5, 1, 2.5)) {
        st <- weibull_survival(t, i, 1)
        plot.dat <- data.frame(t, st)
        p <- ggplot(aes(x = t, y = st), data = plot.dat)
        p <- p + geom_line(size = 1.5)
        p <- p + custom_theme
        p <- p + ylab("S(t)")
        ggsave(paste("Figures/Weibull_surv_p", i, ".svg", sep = ""),
               p,
               device = "svg",
               width = 15.5 / 1.5,
               height = 10.3 / 1.5)
}


## Weibull hazard plots w/ varied p values ##
for (i in c(0.5, 1, 2.5)) {
        st <- weibull_hazard(t, i, 1)
        plot.dat <- data.frame(t, st)
        p <- ggplot(aes(x = t, y = st), data = plot.dat)
        p <- p + geom_line(size = 1.5)
        p <- p + custom_theme
        p <- p + ylab("h(t)")
        ggsave(paste("Figures/Weibull_haz_p", i, ".svg", sep = ""),
               p,
               device = "svg",
               width = 15.5 / 1.5,
               height = 10.3 / 1.5)
}
