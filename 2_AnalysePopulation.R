#!/usr/bin/env Rscript

# This script analyses different populations

# ---- Setup ----

# Cleanup
rm(list = ls())

# Requirements
# devtools::install_bitbucket("hickeyjohnteam/AlphaSimR@devel")
library(package = "AlphaSimR")

# install.packages(pkg = "detectRUNS")
library(package = "detectRUNS")

# install.packages(pkg = "tidyverse")
library(package = "tidyverse")

# install.packages(pkg = "pedigreemm")
library(package = "pedigreemm")

# install.packages("GGally")
library(package = "GGally")

source(file = "Functions.R")

# Arguments
Args = commandArgs(trailingOnly = TRUE)
# Args = c(1, "A")
if (length(Args) < 2) {
  stop("Must provide replicate number as the first argument and scenario as the second argument")
}
Rep = Args[1]
Sce = Args[2]
# Rep = 1
# Sce = "A"

# ---- Data and Constants ----

# Simulated data
ScenarioName = paste("Rep_", Rep, "_Scenario_", Sce, sep = "")
setwd(dir = ScenarioName)
load(file = "Data.RData")

# Pedigree
Ped = read_csv(file = "Pedigree.csv",
               col_types = cols(Gen = col_integer(),
                                IId = col_integer(),
                                FId = col_integer(),
                                MId = col_integer()))

# SNP array map
Map = read_delim(file = "Data.map", col_names = FALSE, delim = " ",
                 col_types = cols(X1 = col_integer(),
                                  X2 = col_character(),
                                  X3 = col_double(),
                                  X4 = col_double())) # have to use double instead of integer due to number formating in file
colnames(Map) = c("Chr", "Loc", "PosGenetic", "Pos")
# head(Map); tail(Map)

# Genome length based on the SNP array
GenomeLengthFromArray = Map %>%
  group_by(Chr) %>%
  summarise(Start  = min(Pos),
            End    = max(Pos),
            Length = End - Start) %>%
  select(Length) %>%
  sum()

# ROH length thresholds based on generations to common ancestor
# Based on the theory of IBD segment length (exponential distribution)
# E(l) = 100 / (2 * GenToCommonAncestor)
GenToCommonAncestor = 1:nGen # c(1:25, 50)
ExpRohLength = 100 / (2 * GenToCommonAncestor)
names(ExpRohLength) = GenToCommonAncestor

# ---- Pedigree (expected) inbreeding ----

MaxGen = max(Ped$Gen)
for (Depth in rev(unique(Ped$Gen))) {
  # Depth = 10
  # Depth =  5
  cat("Pedigree depth", Depth, "generations = Inbreeding from generation", MaxGen - Depth, "\n")
  PedTmp = Ped %>%
    subset(Gen >= (MaxGen - Depth))
  Test = PedTmp$FId %in% PedTmp$IId
  PedTmp$FId[!Test] = NA
  Test = PedTmp$MId %in% PedTmp$IId
  PedTmp$MId[!Test] = NA
  Tmp = pedigree(label = PedTmp$IId,
                 sire  = PedTmp$FId,
                 dam   = PedTmp$MId)
  PedTmp = data_frame(IId = PedTmp$IId,
                      PedInb = inbreeding(ped = Tmp))
  colnames(PedTmp)[2] = paste0("PedInb",
                               formatC(x = Depth, flag = "0", width = nchar(MaxGen)))
  Ped = merge(x = Ped, y = PedTmp, all.x = TRUE)
}
# tail(Ped)

# ---- True (IBD) inbreeding ----

MaxGen = max(Ped$Gen)
for (Depth in rev(unique(Ped$Gen))) {
  # Depth = 5
  # Depth = 0
  cat("Base population removed by", Depth, "generations = Inbreeding from generation", MaxGen - Depth, "\n")
  Tmp = paste0("IbdHaploFrom", formatC(x = MaxGen - Depth, flag = "0", width = nchar(MaxGen)))
  load(file = paste0(Tmp, ".RData"))
  IbdInb = IbdInbreeding(x = get(x = Tmp))
  colnames(IbdInb)[2] = paste0("IbdInb",
                               formatC(x = Depth, flag = "0", width = nchar(MaxGen)))
  Ped = merge(x = Ped, y = IbdInb, all.x = TRUE)
}
# tail(Ped)

# ---- ROH inbreeding ----

ROH <- consecutiveRUNS.run(genotypeFile = "Data.ped", mapFile = "Data.map")
# head(ROH); tail(ROH)
# plot_Runs(runs = ROH)
# plot_StackedRuns(runs = ROH)
MaxGen = max(as.numeric(names(ExpRohLength)))
for (Length in 1L:length(ExpRohLength)) {
  # Length = 1L
  RohInb = RohInbreeding(x = ROH, RohLengthThreshold = ExpRohLength[Length] * 10^6,
                         GenomeLength = GenomeLengthFromArray)
  colnames(RohInb)[1] = "IId"
  RohInb$IId = as.integer(RohInb$IId)
  Gen = as.numeric(names(ExpRohLength)[Length])
  colnames(RohInb)[2] = paste0("RohInb",
                               formatC(x = Gen, flag = "0", width = nchar(MaxGen)))
  Ped = merge(x = Ped, y = RohInb, all.x = TRUE)
}

# ---- Set NA inbreeding to 0 ----

Cols = colnames(Ped)
Cols = Cols[grepl(x = Cols, pattern = "Inb")]
for (Col in Cols) {
  # Col = Cols[1]
  Ped[[Col]] = replace(x = Ped[[Col]], values = 0.0,
                       list = is.na(Ped[[Col]]))
}

# ---- Explore ----

# View(Ped)
if (FALSE) {
  Gen = 10
  ggpairs(Ped[, paste0(c("IbdInb", "RohInb", "PedInb"), Gen)])
}

# TODO: use list vectors? to simplify plotting? or just make the table tall?

cat("DONE\n")

setwd(dir = "..")


# ---- Results data_frame ----

Res = Ped %>%
  subset(Gen == max(Gen))

