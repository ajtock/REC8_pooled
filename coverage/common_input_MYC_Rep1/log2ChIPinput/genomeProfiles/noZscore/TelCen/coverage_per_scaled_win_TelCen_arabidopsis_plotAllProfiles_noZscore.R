#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./coverage_per_scaled_win_TelCen_arabidopsis_plotAllProfiles_noZscore.R 10kb 100 100ths 3 

#winName <- "10kb"
#prop <- 100
#propName <- "100ths" 
#N <- 3 

args <- commandArgs(trailingOnly = T)
winName <- args[1]
prop <- as.numeric(args[2])
propName <- args[3]
N <- as.numeric(args[4])

library(parallel)

profileNames <- c(
                  "log2_ASY1_Rep2_ChIP_REC8_MYC_Rep1_input",
                  "log2_REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input",
                  "log2_MTOPVIB_HA_Rep1_ChIP_REC8_MYC_Rep1_input",
                  "log2_MTOPVIB_HA_Rep2_ChIP_REC8_MYC_Rep1_input",
                  "log2_WT_SPO11oligo_RPI1_WT_nakedDNA_R1"
                 )
profileNamesPlot <- c( 
                      "ASY1",
                      "REC8-HA",
                      "MTOPVIB Rep1",
                      "MTOPVIB Rep2",
                      "SPO11-1-oligos"
                     )
profileColours <- c(
                    "red",
                    "green2",
                    "darkcyan",
                    "cyan",
                    "dodgerblue2"
                   )

makeTransparent <- function(thisColour, alpha = 150)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}
profileColoursTransparent <- sapply(seq_along(profileColours), function(x) {
  makeTransparent(profileColours[x])
})

# Load TelCenMatrix
TelCenDFs <- mclapply(seq_along(profileNames), function(x) {
  read.table(paste0(profileNames[x], "_",
                    winName, "_",
                    propName, "_TelCenMatrix_noZscore.txt"))
}, mc.cores = length(profileNames))

TelCenProfiles <- mclapply(seq_along(TelCenDFs), function(x) {
  as.vector(rowMeans(TelCenDFs[[x]]))
}, mc.cores = length(TelCenDFs))

# Function to plot telomere to centromere (Tel-Cen)
# log2(markChIP/markControl) profiles 
TelCenPlot <- function(xplot,
                       profiles,
                       proportions,
                       proportionsName,
                       profileColours,
                       Ylabel,
                       Ylim,
                       legendLabs,
                       legendLoc) {
  plot(xplot, profiles[[1]], type = "l", lwd = 3, col = profileColours[[1]],
       ylim = Ylim,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "")
  lines(xplot, profiles[[2]], type = "l", lwd = 3, col = profileColours[[2]])
  lines(xplot, profiles[[3]], type = "l", lwd = 3, col = profileColours[[3]])
  lines(xplot, profiles[[4]], type = "l", lwd = 3, col = profileColours[[4]])
  lines(xplot, profiles[[5]], type = "l", lwd = 3, col = profileColours[[5]])
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(1, seq(10, proportions, by = 10)),
       labels = c(expression(italic("TEL")),
                  seq(10, proportions-10, by = 10),
                  expression(italic("CEN"))))
  mtext(side = 1, line = 2.1, cex = 1,
        text = paste0("Scaled windows (", proportionsName, ")"))
  axis(side = 2,
       at = pretty(c(profiles[[1]], profiles[[2]], profiles[[3]], profiles[[4]], profiles[[5]])),
       labels = c("", "-0.4", "", "0.0", "", "0.4", ""),
       cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.1, cex = 1, text = Ylabel, col = "black")
  abline(h = 0, lwd = 1.5, lty = 1)
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c("white"),
         text.col = c(profileColours),
         text.font = c(rep(1, times = 3)),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}


## Calculate moving average of current window,
##### (N/2) previous windows (where N is even) OR
## (N/2)-0.5 previous windows (where N is odd),
## and
##### (N/2) subsequent windows (where N is even) OR
## (N/2)-0.5 subsequent windows (where N is odd)
## (the higher N is, the greater the smoothing)
#stopifnot(N %% 2 != 0) 
#flank <- (N/2)-0.5
## Define MA filter coefficients
#f <- rep(1/N, N)
#
#filt_TelCenProfiles <- mclapply(seq_along(TelCenProfiles), function(x) {
#  filt_TelCenProfile <- stats::filter(x = TelCenProfiles[[x]],
#                                      filter = f,
#                                      sides = 2)
#  filt_TelCenProfile[1:flank] <- filt_TelCenProfile[flank+1]
#  filt_TelCenProfile[(length(filt_TelCenProfile)-flank+1):length(filt_TelCenProfile)] <- filt_TelCenProfile[(length(filt_TelCenProfile)-flank)]
#  filt_TelCenProfile
#}, mc.cores = length(TelCenProfiles))

pdf(paste0("./arabidopsis_log2_ChIPinput_",
           winName, "_",
           propName,
           #propName, "_smooth", N,
           "_TelCenProfile_noZscore_v260919.pdf"),
    height = 3.5, width = 7)
par(mfrow = c(1, 1))
par(mar = c(3.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
TelCenPlot(xplot = 1:length(TelCenProfiles[[1]]),
           profiles = TelCenProfiles,
           proportions = prop,
           proportionsName = propName,
           profileColours = profileColoursTransparent,
           Ylabel = bquote("Log"[2]*"(ChIP/control)"),
           Ylim = c(min(c(TelCenProfiles[[1]], TelCenProfiles[[2]], TelCenProfiles[[3]], TelCenProfiles[[4]], TelCenProfiles[[5]])),
                    max(c(TelCenProfiles[[1]], TelCenProfiles[[2]], TelCenProfiles[[3]], TelCenProfiles[[4]], TelCenProfiles[[5]]))),
           legendLabs = profileNamesPlot,
           legendLoc = "top")
#TelCenPlot(xplot = 1:length(filt_TelCenProfiles[[1]]),
#           profiles = filt_TelCenProfiles,
#           proportions = prop,
#           proportionsName = propName,
#           profileColours = profileColoursTransparent,
#           Ylabel = bquote("Log"[2]*"(ChIP/control)"),
#           Ylim = c(min(c(filt_TelCenProfiles[[1]], filt_TelCenProfiles[[2]], filt_TelCenProfiles[[3]], filt_TelCenProfiles[[4]], filt_TelCenProfiles[[5]])),
#                    max(c(filt_TelCenProfiles[[1]], filt_TelCenProfiles[[2]], filt_TelCenProfiles[[3]], filt_TelCenProfiles[[4]], filt_TelCenProfiles[[5]]))),
#           legendLabs = profileNamesPlot,
#           legendLoc = "top")
dev.off()
