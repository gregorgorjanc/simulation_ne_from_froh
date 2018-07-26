#!/usr/bin/env Rscript

# This script simulates different populations

# ---- Setup ----

# Cleanup
rm(list = ls())

# Requirements
# devtools::install_bitbucket("hickeyjohnteam/AlphaSimR@devel")
library(package = "AlphaSimR")

# install.packages(pkg = "tidyverse")
library(package = "tidyverse")

source(file = "Functions.R")

# Arguments
Args = commandArgs(trailingOnly = TRUE)
# Args = 1
if (length(Args) < 2) {
  stop("Must provide replicate number (1, 2, ...) as the first argument and\n
       scenario (Base, SceA, SceB, SceC, and SceD) as the second argument!")
}
Rep = Args[1]
Sce = Args[2]
# Rep = 1
# Sce = "Base"

# ---- Parameters ----

Species = "CATTLE"
# Species = "GENERIC"

nChr = 30
# nChr =  2
nSnpPerChr = 2000
nQtlPerChr = 1000

nGen       =   25
# nGen       =   10
nIndPerGen = 1000
# nIndPerGen =  100

nSires = 25
nDams = 225
# nDams =  20

# Ne
4 * nSires * nDams / (nSires + nDams)

MeanGInitial = 0
VarGInitial = 1
VarE = 3

# ---- Coalescent and base population ----

# Sce = "Base"
if (Sce == "Base") {
  ScenarioName = paste("Rep_", Rep, sep = "")
  CreateAndSetScenarioFolder(x = ScenarioName)

  FounderPop = runMacs(nInd = nIndPerGen,
                       nChr = nChr,
                       segSites = nSnpPerChr + nQtlPerChr,
                       species = Species)
  SP = SimParam$new(founderPop = FounderPop)
  SP$setGender(gender = "yes_rand")
  SP$addTraitA(nQtlPerChr = nQtlPerChr, mean = MeanGInitial, var = VarGInitial)
  SP$addSnpChip(nSnpPerChr = nSnpPerChr)
  SP$setTrackPed(isTrackPed = TRUE)
  SP$setTrackRec(isTrackRec = TRUE)
  BasePop = newPop(rawPop = FounderPop)
  BasePop = setPheno(pop = BasePop, varE = VarE)

  save.image(file = "Data.RData")

  setwd(dir = "..")
  cat("DONE with", Sce, "\n")
}

# ---- Scenario A: Constant Ne ----

# Sce = "SceA"
if (Sce == "SceA") {
  ScenarioName = paste("Rep_", Rep, "_Scenario_A", sep = "")
  CreateAndSetScenarioFolder(x = ScenarioName)

  load(file = paste("../Rep_", Rep, "/Data.RData", sep = ""))

  Pop = BasePop
  # TODO: save ancestral generations data?
  # Gen = 0
  # writePlink(pop = Pop,
  #            baseName = paste0("Data",
  #                              formatC(x = Gen, flag = "0", width = nchar(nGen))))
  for (Gen in 1:nGen) {
    # Gen = 1
    Sires = selectInd(pop = Pop, nInd = nSires, use = "rand", gender = "M")
    Dams  = selectInd(pop = Pop, nInd = nDams,  use = "rand", gender = "F")
    Pop = randCross2(females = Dams, males = Sires, nCrosses = nIndPerGen,
                     nProgeny = 1, balance = FALSE)
    # TODO: save ancestral generations data?
    # writePlink(pop = Pop,
    #            baseName = paste0("Data",
    #                              formatC(x = Gen, flag = "0", width = nchar(nGen))))
  }
  writePlink(pop = Pop, baseName = "Data")

  Ped = data_frame(Gen = 0,
                   IId = 1:nrow(SP$pedigree),
                   FId = SP$pedigree[, 2],
                   MId = SP$pedigree[, 1])
  Ped$Gen = DetermineGenerationFromBase(x = Ped[, 2:4], Unknown = 0)
  write_csv(x = Ped, path = "Pedigree.csv")

  for (Generation in 0:nGen) {
    # Generation = 0
    PedX = Ped %>%
      mutate(FId = replace(FId, Gen <= Generation, 0),
             MId = replace(MId, Gen <= Generation, 0))
    # View(PedX)
    PedX = PedX %>%
      select(MId, FId) %>%
      as.matrix()
    Tmp = paste0("IbdHaploFrom", formatC(x = Generation, flag = "0", width = nchar(nGen)))
    # TODO: if we remove Pop from the call below we get IBD haplos for the whole ancestral pedigree!
    assign(x = Tmp, value = pullIbdHaplo(pop = Pop, snpChip = 1, pedigree = PedX))
    save(list = Tmp, file = paste0(Tmp, ".RData"))
    rm(list = Tmp); gc()
  }

  save.image(file = "Data.RData")
  setwd(dir = "..")
  cat("DONE with", Sce, "\n")
}

# ---- Scenario B: Constant Ne over time, but 2x of Scenario A due to zero variance in family size ----

# ---- Scenario C: Constant Ne over time, but less than Scenario A due to selection, starts at the same point as Scenario A ----

# ---- Scenario D: Changing Ne over time due to changes in variance of family size, starts at the same point as Scenario A, increases, and then decreases ----

setwd(dir = "..")
