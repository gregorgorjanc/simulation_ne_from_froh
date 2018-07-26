
CreateAndSetScenarioFolder <- function(x) {
  # Create folder and set it as a working directory
  # x - character, folder name
  if (dir.exists(paths = x)) {
    unlink(x = x, recursive = TRUE)
  }
  dir.create(path = x)
  setwd(dir = x)
}

DetermineGenerationFromBase = function(x, Unknown = "0") {
  # Determine how many generations is an individual removed from the base of a pedigree
  # x - pedigree table (tibble or data.frame) with individual, father, and mother columns;
  #     with 0 as unknown parent; the pedigree is assumed sorted such that parents
  #     preceede progeny and that pedigree is extended = that all parents are listed
  #     as individuals)
  # unknown - a value used to denote unknown parent
  if (is.matrix(x)) {
    stop("x must be tibble or data.frame!")
  }
  nInd = nrow(x)
  Generation = rep(x = NA, times = nInd)
  FIDVec = match(x = x[[2L]], table = x[[1L]], nomatch = 0L)
  MIDVec = match(x = x[[3L]], table = x[[1L]], nomatch = 0L)
  FatherIsUnknownVec = x[[2L]] %in% Unknown
  MotherIsUnknownVec = x[[3L]] %in% Unknown
  for (Ind in 1L:nInd) {
    # Ind = 1L
    FID = FIDVec[Ind]
    MID = MIDVec[Ind]
    FatherIsUnknown = FatherIsUnknownVec[Ind]
    MotherIsUnknown = MotherIsUnknownVec[Ind]
    if (FatherIsUnknown & MotherIsUnknown) {
      Generation[Ind] = 0.0
    } else if (!FatherIsUnknown & MotherIsUnknown) {
      Generation[Ind] = (Generation[FID] + 1.0) / 2.0
    } else if (FatherIsUnknown & !MotherIsUnknown) {
      Generation[Ind] = (Generation[MID] + 1.0) / 2.0
    } else {
      Generation[Ind] = ((Generation[FID] + 1.0) + (Generation[MID] + 1.0)) / 2.0
    }
  }
  Generation
}

IbdInbreeding = function(x) {
  # Calculate true (IBD) based inbreeding
  # x - matrix of whole-genome haplotypes (2n)
  Names = rownames(x) %>%
    strsplit(split = "_") %>%
    sapply(FUN = function(z) z[[1L]]) %>%
    as.integer()
  nInd = nrow(x) / 2L
  nLoc = ncol(x)
  Ret = data_frame(IId = rep(x = 0L, times = nInd),
                   IbdInb = 0.0)
  for (Ind in 1L:nInd) {
    # Ind = 1L
    Ret$IId[Ind] = Names[2L * Ind - 2L + 1L]
    # x[2L * Ind - 2L + c(1L, 2L), 1:20]
    Ret$IbdInb[Ind] = sum(x[2L * Ind - 2L + 1L,  ] == x[2L * Ind - 2L + 2L,  ]) / nLoc
  }
  Ret
}

RohInbreeding = function(x, RohLengthThreshold, GenomeLength) {
  # Calculate ROH based inbreeding
  # x - table of ROHs from detectRUNS::consecutiveRUNS.run()
  # RohLengthThreshold - numeric, length of a ROH to be included in analysis
  # GenomeLength - numeric, length of genome to get inbreeding coefficient as percentage
  x %>%
    subset(lengthBps >= RohLengthThreshold) %>%
    group_by(id) %>%
    summarise(RohInb = sum(lengthBps) / GenomeLength)
}