################################################################################
# CROOK MICROPUB DATA ANALYSIS
################################################################################

dat <- read.table("Crook_MP_Data.txt", header = T, sep = "\t", 
                  stringsAsFactors = T) # Read in file

View(dat) # View in spreadsheet-style format

datNoO <- dat[-24, ] # Remove stage 4 outlier

View(datNoO) # View data with outlier removed

################################################################################
# TESTING FOR DIFFERENCES AMONG STAGES
################################################################################

# We will start by fitting an ANOVA model and seeing if the assumptions
# of ANOVA are met.

lMod <- aov(DTCL ~ stage, data = datNoO) # Fit ANOVA model

plot(lMod, which = 1:5) # Check assumptions of ANOVA
# The error variance is not constant
# Residuals don't look horribly non-normal
# No Cook's distances >= 1

shapiro.test(lMod$residuals)
# Reaffirms residuals are more-or-less normal.

# Now we will try log transformation to see if that improves things because
# there are issues with the error variance with the raw data.

logMod <- aov(log(DTCL) ~ stage, data = datNoO)

plot(logMod, which = 1:5)
# Doesn't really improve the error variance.

shapiro.test(logMod$residuals)
# Residuals are more-or-less normal

# For giggles, we will try a square root transformation.

sqrtMod <- aov(sqrt(DTCL) ~ stage, data = datNoO)

plot(sqrtMod, which = 1:5)
# Error variance somewhat better but still doesn't look great.

shapiro.test(sqrtMod$residuals) 
# The residuals are not normal...

# Bottom line, given the small sample size and markedly non-constant error 
# variance we should not run an ANOVA on these data. Let's go ahead and
# make a plot to see if we can get a feel for how similar the shapes of the
# of the distributions are. Admittedly, there's not enough data for each stage
# to do this in a rigorous way.

colo <- gray(0, 0.45) 
# Set color and transparency for plotting characters

par(mgp = c(2.7, 1, 0), cex = 1.35, cex.lab = 1.25, cex.axis = 1.25) 
# Set plotting character, labels, and axis sizes

stripchart(DTCL ~ stage, data = datNoO, vertical = T, pch = 16,
           xlab = "Stage", ylab = "DTC Length (mm)", col = colo,
           method = "jitter")
# Make stripchart

sMeans <- tapply(datNoO$DTCL, INDEX = datNoO$stage, FUN = "mean")
# Calculate means for each stage

segments(x0 = c(0.75, 1.75, 2.75, 3.75), y0 = sMeans, x1 = c(1.25, 2.25, 3.25, 4.25),
         y1 = sMeans, col = "black", lwd = 2.5)
# Overlay stage means on stripchart.

# As stated earlier, there aren't really enough data at each stage to say much 
# about the shapes of the distributions of the various stages. In the interest 
# of not over-thinking the situation we will run a Kruskal-Wallis test.

kwTest <- kruskal.test(DTCL ~ stage, data = datNoO)
kwTest 
# The result is significant. To figure out which stages differ, we will conduct
# pairwise Mann-Whitney U tests and correct for multiple use Holm's sequential
# Bonferroni procedure to control the family-wise error rate (FWER) at the 
# 0.05 level.

kwMC <- pairwise.wilcox.test(x = datNoO$DTCL, g = datNoO$stage, 
                             p.adjust.method = "holm")
kwMC
# This tells us that all of the stages are different from one another, as all
# six adjusted p-values are < 0.05.


################################################################################
# FITTING PIECEWISE MODEL WITH DTCL AS A FUNCTION OF PVL
################################################################################

library(segmented)
# Load the segmented package. If it is not in your library, you will need to 
# install it and all dependencies.

linReg.mod <- lm(DTCL ~ PVL, data = datNoO)
# Fit a simple linear regression model

# Now we need to test for the number of break points, to see how many 
# rate transitions there are.

selgmented(olm = linReg.mod, type = "bic", Kmax = 3)
# The Bayesian Information Criterion (BIC) suggests that a single breakpoint
# is best (note lower BIC indicates better models)

seg.mod <- segmented(obj = linReg.mod)
# Fit the piecewise regression model

# Now we will check the assumptions of the segmented model.
plot(broken.line(seg.mod)$fit, seg.mod$residuals, xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2) 
# there is a bit of a pattern in the residuals suggesting non-constant error
# variance, but given that we're more focused on estimation than hypothesis
# testing here, I don't find it too worrisome.

qqnorm(seg.mod$residuals)
qqline(seg.mod$residuals)
# Residuals don't look all that normal.

shapiro.test(seg.mod$residuals)
# The residuals are not normally distributed, but again, not something that
# is too worrisome given the focus on estimation. It probably does mean that
# the confidence intervals for the breakpoints and slopes should be taken with 
# a healthy grain of salt though, as should p-values associated with 
# t-statistics.

summary(seg.mod)
# Take a look at parameter estimates and fit metrics. The rate transition is
# estimated to occur at a PVL of 0.394 mm. The piecewise model accounts for
# 94.41% of the variation in the data.

slope(seg.mod)
# Take a look at the slopes (rates) for the two line segments. The confidence
# intervals suggest that there is a significant rate change, but these intervals
# should be interpreted with caution (see above).

s1.col <- gray(0, 0.45)
s2.col <- rgb(blue = 1, green = 0, red = 0, alpha = 0.45)
s3.col <- rgb(blue = 0, green = 1, red = 0, alpha = 0.45)
s4.col <- rgb(blue = 0.2, green = 0.6, red = 1, alpha = 0.45)
# Set up plotting colors for the scatter plot.

s.size <- tapply(datNoO$DTCL, datNoO$stage, 'length')
# Determine the sample size for each stage.

s1.c <- rep(s1.col, length.out = s.size[1])
s2.c <- rep(s2.col, length.out = s.size[2])
s3.c <- rep(s3.col, length.out = s.size[3])
s4.c <- rep(s4.col, length.out = s.size[4])
# Continue to set up plotting colors for the scatter plot.

s.c <- c(s1.c, s2.c, s3.c, s4.c)
# Concatenate plotting colors into a single vector for use with 'plot'.

par(mgp = c(2.7, 1, 0), cex = 1.35, cex.lab = 1.25, cex.axis = 1.25)
# Get sizing for the graph right.

plot(datNoO$PVL, datNoO$DTCL, xlab = "Pharynx-vulva length (mm)",
     ylab = "DTC length (mm)", pch = 16, col = s.c)
# Make the scatter plot.

legend("topleft", legend = c("S1", "S2", "S3", "S4"), fill = c(s1.col, s2.col,s3.col, s4.col),
       cex = 0.75)
# Overlay a legend saying what the colors correspond to.

lines.segmented(seg.mod, pch = 16, col = "dark grey", lwd = 1.5)
# Overlay the 95% confidence interval for the breakpoint.

plot.segmented(seg.mod, add = T, col = "black", lwd = 1.5)
# Overlay the piecewise model.

