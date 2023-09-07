

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# source ("mkMarkFile/mk_mark_file_v8.R")
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# refresh workspace
rm (list=ls ())

# load our life-support packages
library (seqinr)

# make it deterministic
set.seed (999)


# ---------------------------------------------------------------------------- #
# load and initialize our data
# ---------------------------------------------------------------------------- #


# set some paths
path.home <- "[path to initial mkMarkFiles directory]"
path.data <- file.path (path.home, "dat")
path.data.neut <- file.path (path.data, "survival")
path.data.map <- file.path (path.data, "map")
path.data.par <- file.path (path.data, "par_scores")
path.data.lineages <- file.path (path.data, "lineages")

# load supporting data:  mark file templates
setwd (path.data.neut)
data.v703 <- read.csv ("v703_survival_wk80_tau_sieve.csv", header=T)
data.v703 <- data.frame (protocol=data.v703[, 1], southAmerica=rep (0, nrow (data.v703)), data.v703[, 2:ncol (data.v703)])
data.v704 <- read.csv ("v704_survival_wk80_tau_sieve.csv", header=T)

# add new columns to the mark file templates for the different sequence types
data.v703 <- data.frame (data.v703[, 1:17], matrix (NA, nrow=nrow (data.v703), ncol=16), data.v703[, 18:ncol (data.v703)])
names (data.v703)[18:33] <- c (gsub ("mf", "ms", names (data.v703)[10:17]), gsub ("mf", "ls", names (data.v703)[10:17]))
data.v704 <- data.frame (data.v704[, 1:17], matrix (NA, nrow=nrow (data.v704), ncol=16), data.v704[, 18:ncol (data.v704)])
names (data.v704)[18:33] <- c (gsub ("mf", "ms", names (data.v704)[10:17]), gsub ("mf", "ls", names (data.v704)[10:17]))

# load supporting data:  HXB2 map
setwd (path.data.map)
hxb2.map.study <- read.csv ("hxb2_study.map", sep="|", header=T)
hxb2.map.lanl <- read.csv ("hxb2_lanl.map", sep="|", header=T)
hxb2.map.catnap <- read.csv ("hxb2_catnap.map", sep="|", header=T)

# load supporting data:  lineages
setwd (path.data.lineages)
lineages <- read.csv ("amp_selected_seqs_v1.csv", header=T)

# load supporting data:  PAR scores
setwd (path.data.par)
par.scores <- read.csv ("par_score_amp_22Nov2021.csv", header=T)

# load supporting data:  MR epitope distances
setwd (file.path (path.data, "epitope_distance"))
epitope.dist.b <- read.csv ("ref_b.tsv", sep="\t", header=T)
epitope.dist.c <- read.csv ("ref_c.tsv", sep="\t", header=T)
epitope.dist.any <- read.csv ("ref_any.tsv", sep="\t", header=T)

# load selected study sequences
setwd (file.path (path.data, "seq/study"))
seq.selected <- read.fasta ("amp_selected_env_aa_v1.fasta", seqtype="AA")

# load supporting data:  reference sequences
setwd (file.path (path.data, "seq/reference") )
refseq.v703 <- read.fasta ("catnap_min_ic80_strains_703_v7.fasta", seqtype="AA")
refseq.v704.lanl <- read.fasta ("lanl_min_ic80_strains_704_v7.fasta", seqtype="AA")
refseq.v704.catnap <- read.fasta ("catnap_min_ic80_strains_704_v7.fasta", seqtype="AA")
refseq.v704 <- c (refseq.v704.catnap, refseq.v704.lanl)
refseq.selections.v703 <- read.csv ("HVTN_703_subset_results_v9.csv", header=T)
refseq.selections.v704 <- read.csv ("HVTN_704_subset_results_v9.csv", header=T)[, -6]

# load supporting data:  site information
setwd (file.path (path.data, "lookup"))
site.info <- read.csv ("pub_id_site_lookup_cam.csv", header=T)

# load supporting data:  miscellany
setwd (file.path (path.data, "misc"))
env.map <- read.csv ("env_locations.csv", header=T)
sites <- read.csv ("sites_residues.csv", header=T)
cd4.vrc01.bsites <- read.csv ("bsites_cd4_vrc01.dat", header=F)[, 1]
geo <- read.csv ("pub_id_site_lookup_cam_stratasafrica.csv", header=T)
subtypes.table <- read.csv ("subtype_final.csv", header=T)

# make a filter to keep pre-selected sites that aren't hypervariable
sites.filter.nonhypervar <- !(sites[, 1] %in% c (133:154, 185, 187:192, 393:416, 460:466))

# QC:  confirm we have a PAR score for each selected sequence
#table (as.vector (na.omit (lineages$seqid.aa)) %in% par.scores$seqname)

# define some other variables to collect
sites.pngs <- c (156, 229, 234, 616, 824)
coords.gp120 <- env.map[env.map[, 1] == "GP120", 2:3]
coords.v1v2 <- c (env.map[env.map[, 1] == "V1_loop", 2], env.map[env.map[, 1] == "V2_loop", 3])
coords.v5 <- env.map[env.map[, 1] == "V5_loop", 2:3]

# which sites are we including in:  "The set of 17 AA positions found to be 
# "important either for predicting neutralization resistance or for PNGS, 
# i.e., the union of positions listed in Table 2 and those listed as 
# PNGS-related sites in Section 5.2.1"
important.sites.all <- sort (unique (c (sites.pngs, sites[sites.filter.nonhypervar, 1])))

# load and initialize our z-space table
setwd (file.path (path.data, "zspace"))
zspace.table <- read.csv ("wold_z5_space_diff_v3.csv", header=T)
rownames (zspace.table) <- zspace.table[, 1]
zspace.table <- zspace.table[, -1]
colnames (zspace.table) <- rownames (zspace.table)

# make a list of pubids
pubids.pooled <- sort (unique (substr (names (seq.selected), 2, 9)))
pubids.pooled <- gsub ("_", "-", pubids.pooled)
pubids.v703 <- pubids.pooled[substr (pubids.pooled, 1, 3) == "703"]
pubids.v704 <- pubids.pooled[substr (pubids.pooled, 1, 3) == "704"]


# ---------------------------------------------------------------------------- #
# FXN Depot
# ---------------------------------------------------------------------------- #


# create a function that counts the number of PNGSes in a sequence, using the
# standard PNGS motif:  [N][!P][S|T]
num.pngs <- function (seq) {
  count <- 0
  for (start in 1:(length (seq) - 2)) {
    count <- count + as.numeric (seq[start] == "N" & seq[start + 1] != "P"& seq[start + 2] %in% c ("S", "T"))
  }
  return (count)
}

# a function for calculating z-space-weighted Hamming distances between two 
# sequence vectors
hdist.zspace <- function (seq.tmp, insert.tmp) {
  dist.tmp <- 0
  seq.tmp[is.na (seq.tmp)] <- "-"
  insert.tmp[is.na (insert.tmp)] <- "-"
  for (pos.tmp in 1:length (insert.tmp)) {
    dist.tmp <- dist.tmp + zspace.table[seq.tmp[pos.tmp], insert.tmp[pos.tmp]]
  }
  return (dist.tmp)
}

# confirm whether or not we pass the minvar threshold
is.minvar <- function (residue.vector, threshold=4) {
  table.tmp <- table (residue.vector)
  return (max (table.tmp) <= length (residue.vector) - threshold)
}

# return the logit of a probability
logit <- function (p) {
  return (log (p / (1 - p)))
}

# parscore3:  ordered category of predicted IC80 ("<1", 1-3 ("[1,3]"), or ">3" ug/ml)
parscore3 <- function (par) {
  if (par < 1) {
    return ("<1")
  } else if (par >= 1 & par <= 3) {
    return ("[1,3]")
  } else if (par > 3) {
    return (">3")
  }
}

# hiv1eventparscore3:  indicator of HIV infected primary case by parscore3
hiv1eventparscore3 <- function (par) {
  if (par < 1) {
    return (1)
  } else if (par >= 1 & par <= 3) {
    return (2)
  } else if (par > 3) {
    return (3)
  }
}

# is a given 3mer a sequon?
is.sequon <- function (triplet) {
  return (as.numeric ((triplet[1] == "N" & triplet[2] != "P" & triplet[3] %in% c ("S", "T"))))
}

# get the 3mer starting at a given position from a given sequence
get.sequon <- function (inseq, start.pos, hxb2.map) {
  index.tmp <- hxb2.map[hxb2.map[, 2] == start.pos, 1]
  subseq <- inseq[index.tmp:length (inseq)]
  return (subseq[subseq != "-"][1:3])
}

# wrapper for get.sequon to run on multiple sites
get.all.sequons <- function (inseq, sites, hxb2.map) {
  results <- NULL
  for (site in sites) {
    results <- append (results, get.sequon (inseq, site, hxb2.map))
  }
  return (results)
}

# get.subseq:  given a sequence as a vector, and start/stop positions, return
#              the subsequence corresponding to those positions
get.subseq <- function (seq, start, stop) {
  subseq <- seq[start:stop]
  return (subseq[subseq != "-"])
}


# ---------------------------------------------------------------------------- #
# SIDE WORK:  for each participant, select our most-frequent, most-sensitive 
#             and least-sensitive sequence
# ---------------------------------------------------------------------------- #


# initialize our results object
rep.seqs <- data.frame ()

# loop through each pubid
for (pubid in sort (unique (lineages$pubid))) {

  # subset our sequences
  seq.tmp <- lineages[lineages$pubid == pubid, ]

  # add our associated PAR scores
  seq.tmp <- merge (seq.tmp, par.scores[par.scores$seqname %in% seq.tmp$seqid.aa, c (4, 6)], by.x="seqid.aa", by.y="seqname")

  # select the sequences for our three categories:  most-frequent
  mf <- seq.tmp[seq.tmp$seq.freq.aa == max (seq.tmp$seq.freq.aa), ]
  mf <- mf[mf$hdist.nt.from.maxseq == min (mf$hdist.nt.from.maxseq), "seqid.aa"]

  # select the sequences for our three categories:  most-sensitive
  ms <- seq.tmp[seq.tmp$pred.ic80 == min (seq.tmp$pred.ic80), ]
  ms <- sample (ms[ms$seq.freq.aa == max (ms$seq.freq.aa), "seqid.aa"], 1)

  # select the sequences for our three categories:  least-sensitive
  ls <- seq.tmp[seq.tmp$pred.ic80 == max (seq.tmp$pred.ic80), ]
  ls <- sample (ls[ls$seq.freq.aa == max (ls$seq.freq.aa), "seqid.aa"], 1)

  # compile our results
  rep.seqs <- rbind (rep.seqs, data.frame (pubid, mf, ms, ls))
}

# add column for study information
rep.seqs$study <- substr (rep.seqs[, 1], 1, 4)
rep.seqs <- rep.seqs[, c (5, 1:4)]


# ---------------------------------------------------------------------------- #
# SIDE WORK:  create a table of which of the fifty sites passed the minimum 
#             variability criterion for each selected sequence, for each study
# ---------------------------------------------------------------------------- #


# create a dictionary mapping our full sequence names with the shorter 
# sequence IDs
seq.id.lookup <- list ()
for (seqname in names (seq.selected)) {
  bits <- unlist (strsplit (seqname, split="|", fixed=T))
  seq.id.lookup[[bits[1]]] <- seqname
}

# initialize our results object
keepers <- list ()
keepers[["v703"]] <- list ()
keepers[["v704"]] <- list ()
keepers[["pooled"]] <- list ()

# loop through each site
for (site in cd4.vrc01.bsites) {

  # get the alignment position for this site
  index <- hxb2.map.study[hxb2.map.study[, 2] == site, "posNum"]

  # loop through each sequence type
  for (seq.type in c ("mf", "ms", "ls")) {

    # initialize some tally objects
    residues.v703 <- NULL
    residues.v704 <- NULL
    residues.pooled <- NULL

    # loop through each participant:  HVTN 703
    for (pubid in rep.seqs[rep.seqs$study == "V703", "pubid"]) {
    
      # identify the residue for each sequence and add them to the list
      seqname.tmp <- rep.seqs[rep.seqs$pubid == pubid, seq.type]
      residues.v703 <- append (residues.v703, seq.selected[[seq.id.lookup[[seqname.tmp]]]][index])
      residues.pooled <- append (residues.pooled, seq.selected[[seq.id.lookup[[seqname.tmp]]]][index])
    }

    # loop through each participant:  HVTN 704
    for (pubid in rep.seqs[rep.seqs$study == "V704", "pubid"]) {
    
      # identify the residue for each sequence and add them to the list
      seqname.tmp <- rep.seqs[rep.seqs$pubid == pubid, seq.type]
      residues.v704 <- append (residues.v704, seq.selected[[seq.id.lookup[[seqname.tmp]]]][index])
      residues.pooled <- append (residues.pooled, seq.selected[[seq.id.lookup[[seqname.tmp]]]][index])
    }
    
    # initialize our results objects, if necessary
    if (!seq.type %in% names (keepers[["v703"]])) {
      keepers[["v703"]][[seq.type]] <- is.minvar (residues.v703)
      keepers[["v704"]][[seq.type]] <- is.minvar (residues.v704)
      keepers[["pooled"]][[seq.type]] <- is.minvar (residues.pooled)
    } else {
      keepers[["v703"]][[seq.type]] <- append (keepers[["v703"]][[seq.type]], is.minvar (residues.v703))
      keepers[["v704"]][[seq.type]] <- append (keepers[["v704"]][[seq.type]], is.minvar (residues.v704))
      keepers[["pooled"]][[seq.type]] <- append (keepers[["pooled"]][[seq.type]], is.minvar (residues.pooled))
    }
  }
}

# compile results into dataframe
sites.included <- data.frame (position=cd4.vrc01.bsites,
                              v703.mf=keepers[["v703"]][["mf"]],
                              v703.ms=keepers[["v703"]][["ms"]],
                              v703.ls=keepers[["v703"]][["ls"]],
                              v704.mf=keepers[["v704"]][["mf"]],
                              v704.ms=keepers[["v704"]][["ms"]],
                              v704.ls=keepers[["v704"]][["ls"]],
                              pooled.mf=keepers[["pooled"]][["mf"]],
                              pooled.ms=keepers[["pooled"]][["ms"]],
                              pooled.ls=keepers[["pooled"]][["ls"]])


# ---------------------------------------------------------------------------- #
# SIDE WORK:  determine which residues are seen at each VRC01/CD4 binding site
# ---------------------------------------------------------------------------- #


# initialize a list for our results
sites.residues <- list ()

# loop through each site
for (site.tmp in cd4.vrc01.bsites) {

  # sort out our alignment index
  index.study.tmp <- hxb2.map.study[hxb2.map.study[, 2] == site.tmp, 1]

  # compile all residues observed at this site
  residues.tmp <- NULL
  for (seqnum.tmp in 1:length (seq.selected)) {
    residues.tmp <- append (residues.tmp, seq.selected[[seqnum.tmp]][index.study.tmp])
  }

  # sort out our individual observed residues at this site
  sites.residues[[site.tmp]] <- sort (unique (residues.tmp))
}


# ---------------------------------------------------------------------------- #
# CONFIG:  which study are we doing this for?
# ---------------------------------------------------------------------------- #


# set up our study-specific data objects accordingly
data <- rbind (data.v703, data.v704)
refseq.selections <- rbind (refseq.selections.v703, refseq.selections.v704)
refseqs <- append (refseq.v703, refseq.v704)
pubids <- pubids.pooled


# ---------------------------------------------------------------------------- #
# take care of business
#
# data required:
#   -- Tier 1:  (1a) logit predicted probability IC80 > 1 ug/ml
#               (1b) log base 10 predicted IC80
#               (1c) Ordered category of predicted IC80 < 1 vs. 1-3 vs. > 3 ug/ml
#               (2)  Levels of all important features predicting neutralization 
#                    sensitivity
#   -- Tier 2:  (3)  AA site scanning: VRC01+CD4bs:  features of all 
#                    sufficiently-variable positions in the VRC01 footprint or 
#                    the CD4 binding region
#               (4)  Distances to most sensitive subtype-matched sequence:
#                    (a) Hamming dists all important features in (2)
#                    (b) Hamming dists in subset of features in (2) with a sieve effect (to be done later)
#                    (c) Hamming dists in all VRC01+CD4bs in (3)
#                    (d) Hamming dists in subset of VRC01+CD4bs in (3) with a sieve effect (to be done later)
#
# ---------------------------------------------------------------------------- #


# format our data frames of sequence residue data
residues.preselect <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=nrow (sites) * 3))
names (residues.preselect) <- paste (as.vector (sapply (paste ("hxb2", sites[, 1], sites[, 2], sep="."), rep, 3)), c ("mf", "ms", "ls"), sep=".")
residues.bsites <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=length (cd4.vrc01.bsites) * 3))
names (residues.bsites) <- paste (as.vector (sapply (paste ("hxb2", cd4.vrc01.bsites, "ref", sep="."), rep, 3)), c ("mf", "ms", "ls"), sep=".")
pngs <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=(length (sites.pngs) * 3)))
names (pngs) <- paste (as.vector (sapply (paste ("hxb2", sites.pngs, "pngs", sep="."), rep, 3)), c ("mf", "ms", "ls"), sep=".")
virgeom <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=(7 * 3)))
virgeom.names <- c ("length.gp120", "length.v1v2", "length.v5", "num.pngs.gp120", "num.pngs.v1v2", "num.pngs.v5", "num.cysteine.gp120")
num.founders <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=2))
names (num.founders) <- c ("num.founders.all", "num.founders.tfl")
names (virgeom) <- paste (as.vector (sapply (virgeom.names, rep, 3)), c ("mf", "ms", "ls"), sep=".")
hdist <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=(4 * 3)))
names (hdist) <- c ("hdist.zspace.sites.preselect.all.mf", "hdist.zspace.sites.preselect.all.ms", "hdist.zspace.sites.preselect.all.ls",
                    "hdist.zspace.sites.preselect.sieve.mf", "hdist.zspace.sites.preselect.sieve.ms", "hdist.zspace.sites.preselect.sieve.ls",
                    "hdist.zspace.sites.binding.all.mf", "hdist.zspace.sites.binding.all.ms", "hdist.zspace.sites.binding.all.ls",
                    "hdist.zspace.sites.binding.sieve.mf", "hdist.zspace.sites.binding.sieve.ms", "hdist.zspace.sites.binding.sieve.ls")
epitope.dist <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=(3 * 3)))
names (epitope.dist) <- c ("epitope.dist.b.mf", "epitope.dist.b.ms", "epitope.dist.b.ls",
                           "epitope.dist.c.mf", "epitope.dist.c.ms", "epitope.dist.c.ls",
                           "epitope.dist.any.mf", "epitope.dist.any.ms", "epitope.dist.any.ls")
founders.sensitive.indicators <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=4))
names (founders.sensitive.indicators) <- c ("all.founders.sensitive", "at.least.one.founder.resistant", "all.founders.resistant", "all.founders.all.sens.resist")
glycosite.230 <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=6))
names (glycosite.230) <- c ("hxb2.230.N.mf", "hxb2.230.N.ms", "hxb2.230.N.ls", "hxb2.230.pngs.mf", "hxb2.230.pngs.ms", "hxb2.230.pngs.ls")
subtype <- rep (NA, nrow (data))
country <- rep (NA, nrow (data))

# update our default hiv1eventparscore3 values
data$hiv1eventparscore3.mf[data$hiv1event == 0] <- 0
data$hiv1eventparscore3.ms[data$hiv1event == 0] <- 0
data$hiv1eventparscore3.ls[data$hiv1event == 0] <- 0

# loop through each participant
for (pubid.num in 1:length (pubids)) {

  # set up some convenience variables
  pubid.raw.tmp <- pubids[pubid.num]
  pubid.seq.tmp <- paste0 ("V", gsub ("-", "_", pubid.raw.tmp))
  lineages.tmp <- lineages[lineages$pubid == pubid.seq.tmp, ]
  par.tmp <- par.scores[par.scores$seqname %in% lineages.tmp$seqid.aa, ]
  lineages.tmp <- merge (lineages.tmp, par.tmp[, c ("seqname", "pred.ic80")], by.x="seqid.aa", by.y="seqname")

  # determine the row of the results we want to modify
  results.row <- data$pub_id == pubid.raw.tmp

  # identify our representative sequences
  seqid.mf <- sample (lineages.tmp[lineages.tmp$seq.freq.aa == max (lineages.tmp$seq.freq.aa), "seqid.aa"], 1)
  seqid.ms <- sample (lineages.tmp[lineages.tmp$pred.ic80 == min (lineages.tmp$pred.ic80), "seqid.aa"], 1)
  seqid.ls <- sample (lineages.tmp[lineages.tmp$pred.ic80 == max (lineages.tmp$pred.ic80), "seqid.aa"], 1)

  # note the actual representative sequences
  seq.mf <- as.vector (seq.selected[[seq.id.lookup[[seqid.mf]]]])
  seq.ms <- as.vector (seq.selected[[seq.id.lookup[[seqid.ms]]]])
  seq.ls <- as.vector (seq.selected[[seq.id.lookup[[seqid.ls]]]])

  # get the PAR scores for our selected sequences
  par.scores.mf <- par.scores[par.scores$seqname == seqid.mf, 5:6]
  par.scores.ms <- par.scores[par.scores$seqname == seqid.ms, 5:6]
  par.scores.ls <- par.scores[par.scores$seqname == seqid.ls, 5:6]

  # identify our reference sequence and the alignment map for the source
  seqid.ref <- refseq.selections[refseq.selections$pubid == pubid.raw.tmp, "LANL.Min.Seqname"]
  seq.ref <- as.vector (refseqs[[seqid.ref]])
  seq.ref.source <- refseq.selections[refseq.selections$pubid == pubid.raw.tmp, "Reference.Source"]
  if (seq.ref.source == "CATNAP") {
    hxb2.map.ref <- hxb2.map.catnap
  } else if (seq.ref.source == "LANL") {
    hxb2.map.ref <- hxb2.map.lanl
  }

  # TIER 1 (1a/b/c):  most-frequent
  data[results.row, "seqname.mf"] <- seqid.mf
  data[results.row, "pred.ic80.mf"] <- par.scores.mf["pred.ic80"]
  data[results.row, "pred.prob.sens.mf"] <- par.scores.mf["pred.prob.sens"]
  data[results.row, "pred.prob.res.mf"] <- 1 - par.scores.mf["pred.prob.sens"]
  data[results.row, "parscore1.mf"] <- logit (1 - par.scores.mf["pred.prob.sens"])
  data[results.row, "parscore2.mf"] <- log10 (par.scores.mf["pred.ic80"])
  data[results.row, "parscore3.mf"] <- parscore3 (par.scores.mf["pred.ic80"])
  data[results.row, "hiv1eventparscore3.mf"] <- hiv1eventparscore3 (par.scores.mf["pred.ic80"])

  # TIER 1 (1a/b/c):  most-sensitive
  data[results.row, "seqname.ms"] <- seqid.ms
  data[results.row, "pred.ic80.ms"] <- par.scores.ms["pred.ic80"]
  data[results.row, "pred.prob.sens.ms"] <- par.scores.ms["pred.prob.sens"]
  data[results.row, "pred.prob.res.ms"] <- 1 - par.scores.ms["pred.prob.sens"]
  data[results.row, "parscore1.ms"] <- logit (1 - par.scores.ms["pred.prob.sens"])
  data[results.row, "parscore2.ms"] <- log10 (par.scores.ms["pred.ic80"])
  data[results.row, "parscore3.ms"] <- parscore3 (par.scores.ms["pred.ic80"])
  data[results.row, "hiv1eventparscore3.ms"] <- hiv1eventparscore3 (par.scores.ms["pred.ic80"])

  # TIER 1 (1a/b/c):  least-sensitive
  data[results.row, "seqname.ls"] <- seqid.ls
  data[results.row, "pred.ic80.ls"] <- par.scores.ls["pred.ic80"]
  data[results.row, "pred.prob.sens.ls"] <- par.scores.ls["pred.prob.sens"]
  data[results.row, "pred.prob.res.ls"] <- 1 - par.scores.ls["pred.prob.sens"]
  data[results.row, "parscore1.ls"] <- logit (1 - par.scores.ls["pred.prob.sens"])
  data[results.row, "parscore2.ls"] <- log10 (par.scores.ls["pred.ic80"])
  data[results.row, "parscore3.ls"] <- parscore3 (par.scores.ls["pred.ic80"])
  data[results.row, "hiv1eventparscore3.ls"] <- hiv1eventparscore3 (par.scores.ls["pred.ic80"])

  # TIER 1 (2):  residue matches (to resistant residues) at pre-specified sites
  for (site.num in 1:nrow (sites)) {
    site.tmp <- sites[site.num, 1]
    index.tmp <- hxb2.map.study[hxb2.map.study[, 2] == site.tmp, 1]
    resist.residue.tmp <- sites[site.num, 2]
    seq.residue.mf.tmp <- seq.mf[index.tmp]
    seq.residue.ms.tmp <- seq.ms[index.tmp]
    seq.residue.ls.tmp <- seq.ls[index.tmp]
    seq.matches.tmp <- c (as.numeric (seq.residue.mf.tmp == resist.residue.tmp),
                          as.numeric (seq.residue.ms.tmp == resist.residue.tmp),
                          as.numeric (seq.residue.ls.tmp == resist.residue.tmp))
    residues.preselect[results.row, (((site.num - 1) * 3) + 1):((((site.num - 1) * 3) + 1) + 2)] <- seq.matches.tmp
  }

  # TIER 1 (2):  potential PNGSes at pre-specified sites
  for (site.num in 1:length (sites.pngs)) {
    site.tmp <- sites.pngs[site.num]
    index.tmp <- hxb2.map.study[hxb2.map.study[, 2] == site.tmp, 1]
    subseq <- seq.mf[index.tmp:length (seq.mf)]
    subseq <- subseq[subseq != "-"][1:3]
    is.sequon.mf <- is.sequon (subseq)
    subseq <- seq.ms[index.tmp:length (seq.ms)]
    subseq <- subseq[subseq != "-"][1:3]
    is.sequon.ms <- is.sequon (subseq)
    subseq <- seq.ls[index.tmp:length (seq.ls)]
    subseq <- subseq[subseq != "-"][1:3]
    is.sequon.ls <- is.sequon (subseq)
    pngs[results.row, (((site.num - 1) * 3) + 1):((((site.num - 1) * 3) + 1) + 2)] <- c (is.sequon.mf, is.sequon.ms, is.sequon.ls)
  }

  # TIER 1 (2) determine the preselected viral geometry features:  gp120
  subseq.gp120 <- get.subseq (seq.mf, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.gp120[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.gp120[2]), 1])
  subseq.v1v2 <- get.subseq (seq.mf, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v1v2[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v1v2[2]), 1])
  subseq.v5 <- get.subseq (seq.mf, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v5[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v5[2]), 1])
  virgeom[results.row, "length.gp120.mf"] <- length (subseq.gp120)
  virgeom[results.row, "length.v1v2.mf"] <- length (subseq.v1v2)
  virgeom[results.row, "length.v5.mf"] <- length (subseq.v5)
  virgeom[results.row, "num.pngs.gp120.mf"] <- num.pngs (subseq.gp120)
  virgeom[results.row, "num.pngs.v1v2.mf"] <- num.pngs (subseq.v1v2)
  virgeom[results.row, "num.pngs.v5.mf"] <- num.pngs (subseq.v5)
  virgeom[results.row, "num.cysteine.gp120.mf"] <- sum (subseq.gp120 == "C")
  subseq.gp120 <- get.subseq (seq.ms, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.gp120[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.gp120[2]), 1])
  subseq.v1v2 <- get.subseq (seq.ms, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v1v2[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v1v2[2]), 1])
  subseq.v5 <- get.subseq (seq.ms, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v5[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v5[2]), 1])
  virgeom[results.row, "length.gp120.ms"] <- length (subseq.gp120)
  virgeom[results.row, "length.v1v2.ms"] <- length (subseq.v1v2)
  virgeom[results.row, "length.v5.ms"] <- length (subseq.v5)
  virgeom[results.row, "num.pngs.gp120.ms"] <- num.pngs (subseq.gp120)
  virgeom[results.row, "num.pngs.v1v2.ms"] <- num.pngs (subseq.v1v2)
  virgeom[results.row, "num.pngs.v5.ms"] <- num.pngs (subseq.v5)
  virgeom[results.row, "num.cysteine.gp120.ms"] <- sum (subseq.gp120 == "C")
  subseq.gp120 <- get.subseq (seq.ls, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.gp120[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.gp120[2]), 1])
  subseq.v1v2 <- get.subseq (seq.ls, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v1v2[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v1v2[2]), 1])
  subseq.v5 <- get.subseq (seq.ls, hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v5[1]), 1], hxb2.map.study[hxb2.map.study[, 2] == as.character (coords.v5[2]), 1])
  virgeom[results.row, "length.gp120.ls"] <- length (subseq.gp120)
  virgeom[results.row, "length.v1v2.ls"] <- length (subseq.v1v2)
  virgeom[results.row, "length.v5.ls"] <- length (subseq.v5)
  virgeom[results.row, "num.pngs.gp120.ls"] <- num.pngs (subseq.gp120)
  virgeom[results.row, "num.pngs.v1v2.ls"] <- num.pngs (subseq.v1v2)
  virgeom[results.row, "num.pngs.v5.ls"] <- num.pngs (subseq.v5)
  virgeom[results.row, "num.cysteine.gp120.ls"] <- sum (subseq.gp120 == "C")

  # TIER 2 (3):  residue matches (to reference sequences) at CD4/VRC01 binding sites
  for (site.num in 1:length (cd4.vrc01.bsites)) {

    # identify the site
    site.tmp <- cd4.vrc01.bsites[site.num]

    # sort out our alignment indices
    index.study.tmp <- hxb2.map.study[hxb2.map.study[, 2] == site.tmp, 1]
    index.catnap.tmp <- hxb2.map.catnap[hxb2.map.catnap[, 2] == site.tmp, 1]
    index.lanl.tmp <- hxb2.map.lanl[hxb2.map.lanl[, 2] == site.tmp, 1]
    if (seq.ref.source == "CATNAP") {
      index.source.tmp <- index.catnap.tmp
    } else if (seq.ref.source == "LANL") {
      index.source.tmp <- index.lanl.tmp
    }

    # capture our residues
    seq.residue.mf.tmp <- seq.mf[index.study.tmp]
    seq.residue.ms.tmp <- seq.ms[index.study.tmp]
    seq.residue.ls.tmp <- seq.ls[index.study.tmp]
    seq.residue.ref.tmp <- seq.ref[index.source.tmp]

    # calculate and immortalize our matches
    seq.matches.tmp <- c (as.numeric (seq.residue.mf.tmp == seq.residue.ref.tmp),
                          as.numeric (seq.residue.ms.tmp == seq.residue.ref.tmp),
                          as.numeric (seq.residue.ls.tmp == seq.residue.ref.tmp))
    residues.bsites[results.row, (((site.num - 1) * 3) + 1):((((site.num - 1) * 3) + 1) + 2)] <- seq.matches.tmp
  }

  # TIER 2 (4a):  weighted Hamming distances of all pre-selected sites
  seq.mf.sites <- c (seq.mf[hxb2.map.study[, 2] %in% sites[sites.filter.nonhypervar, 1]], get.all.sequons (seq.mf, sites.pngs, hxb2.map.study))
  seq.ms.sites <- c (seq.ms[hxb2.map.study[, 2] %in% sites[sites.filter.nonhypervar, 1]], get.all.sequons (seq.ms, sites.pngs, hxb2.map.study))
  seq.ls.sites <- c (seq.ls[hxb2.map.study[, 2] %in% sites[sites.filter.nonhypervar, 1]], get.all.sequons (seq.ls, sites.pngs, hxb2.map.study))
  seq.ref.sites <- c (seq.ref[hxb2.map.ref[, 2] %in% sites[sites.filter.nonhypervar, 1]], get.all.sequons (seq.ref, sites.pngs, hxb2.map.ref))
  hdist[results.row, "hdist.zspace.sites.preselect.all.mf"] <- hdist.zspace (seq.mf.sites, seq.ref.sites)
  hdist[results.row, "hdist.zspace.sites.preselect.all.ms"] <- hdist.zspace (seq.ms.sites, seq.ref.sites)
  hdist[results.row, "hdist.zspace.sites.preselect.all.ls"] <- hdist.zspace (seq.ls.sites, seq.ref.sites)

  # TIER 2 (4b):  Hamming dists in subset of features in (2) with a sieve effect
  #   TO DO LATER

  # TIER 2 (4c):  weighted Hamming distances of all CD4/VRC01 binding sites
  seq.mf.bsites <- seq.mf[hxb2.map.study[, 2] %in% cd4.vrc01.bsites]
  seq.ms.bsites <- seq.ms[hxb2.map.study[, 2] %in% cd4.vrc01.bsites]
  seq.ls.bsites <- seq.ls[hxb2.map.study[, 2] %in% cd4.vrc01.bsites]
  seq.ref.bsites <- seq.ref[hxb2.map.ref[, 2] %in% cd4.vrc01.bsites]
  hdist[results.row, "hdist.zspace.sites.binding.all.mf"] <- hdist.zspace (seq.mf.bsites, seq.ref.bsites)
  hdist[results.row, "hdist.zspace.sites.binding.all.ms"] <- hdist.zspace (seq.ms.bsites, seq.ref.bsites)
  hdist[results.row, "hdist.zspace.sites.binding.all.ls"] <- hdist.zspace (seq.ls.bsites, seq.ref.bsites)

  # TIER 2 (4d):  Hamming dists in subset of VRC01+CD4bs in (3) with a sieve effect
  #   TO DO LATER

  # EXTRA:  indicators of founder neutralization sensitivity, i.e.:
  #         -- indicator of all founders (i.e., all founder cluster-specific 
  #            mindist sequences) being predicted sensitive,
  #         -- indicator of at least one founder being predicted resistant,
  #         -- indicator of all founders being predicted resistant.
  founders.sensitive.indicators[results.row, "all.founders.sensitive"] <- as.numeric (data[results.row, "pred.ic80.mf"] < 1 & data[results.row, "pred.ic80.ms"] < 1 & data[results.row, "pred.ic80.ls"] < 1)
  founders.sensitive.indicators[results.row, "at.least.one.founder.resistant"] <- as.numeric (sum (c (data[results.row, "pred.ic80.mf"] >= 1, data[results.row, "pred.ic80.ms"] >= 1, data[results.row, "pred.ic80.ls"] >= 1)) >= 1)
  founders.sensitive.indicators[results.row, "all.founders.resistant"] <- as.numeric (sum (c (data[results.row, "pred.ic80.mf"] >= 1, data[results.row, "pred.ic80.ms"] >= 1, data[results.row, "pred.ic80.ls"] >= 1)) == 3)
  if (founders.sensitive.indicators[results.row, "all.founders.sensitive"] == 1) {
    founders.sensitive.indicators[results.row, "all.founders.all.sens.resist"] <- 1
  } else if (founders.sensitive.indicators[results.row, "all.founders.resistant"] == 1) {
    founders.sensitive.indicators[results.row, "all.founders.all.sens.resist"] <- 0
  }

  # EXTRA:  details about glycosite 230
  glycosite.230[results.row, "hxb2.230.N.mf"] <- as.numeric (seq.mf[hxb2.map.study[hxb2.map.study[, 2] == "230", 1]] == "N")
  glycosite.230[results.row, "hxb2.230.N.ms"] <- as.numeric (seq.ms[hxb2.map.study[hxb2.map.study[, 2] == "230", 1]] == "N")
  glycosite.230[results.row, "hxb2.230.N.ls"] <- as.numeric (seq.ls[hxb2.map.study[hxb2.map.study[, 2] == "230", 1]] == "N")
  glycosite.230[results.row, "hxb2.230.pngs.mf"] <- is.sequon (get.sequon (seq.mf, "230", hxb2.map.study))
  glycosite.230[results.row, "hxb2.230.pngs.ms"] <- is.sequon (get.sequon (seq.ms, "230", hxb2.map.study))
  glycosite.230[results.row, "hxb2.230.pngs.ls"] <- is.sequon (get.sequon (seq.ls, "230", hxb2.map.study))

  # LATE COMERS:  epitope distances provided by Morgane Rolland
  if (length (epitope.dist.b[epitope.dist.b$strains == seqid.mf, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.b.mf"] <- epitope.dist.b[epitope.dist.b$strains == seqid.mf, "ep_distance"]
  }
  if (length (epitope.dist.b[epitope.dist.b$strains == seqid.ms, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.b.ms"] <- epitope.dist.b[epitope.dist.b$strains == seqid.ms, "ep_distance"]
  }
  if (length (epitope.dist.b[epitope.dist.b$strains == seqid.ls, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.b.ls"] <- epitope.dist.b[epitope.dist.b$strains == seqid.ls, "ep_distance"]
  }
  if (length (epitope.dist.c[epitope.dist.c$strains == seqid.mf, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.c.mf"] <- epitope.dist.c[epitope.dist.c$strains == seqid.mf, "ep_distance"]
  }
  if (length (epitope.dist.c[epitope.dist.c$strains == seqid.ms, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.c.ms"] <- epitope.dist.c[epitope.dist.c$strains == seqid.ms, "ep_distance"]
  }
  if (length (epitope.dist.c[epitope.dist.c$strains == seqid.ls, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.c.ls"] <- epitope.dist.c[epitope.dist.c$strains == seqid.ls, "ep_distance"]
  }
  if (length (epitope.dist.any[epitope.dist.any$strains == seqid.mf, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.any.mf"] <- epitope.dist.any[epitope.dist.any$strains == seqid.mf, "ep_distance"]
  }
  if (length (epitope.dist.any[epitope.dist.any$strains == seqid.ms, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.any.ms"] <- epitope.dist.any[epitope.dist.any$strains == seqid.ms, "ep_distance"]
  }
  if (length (epitope.dist.any[epitope.dist.any$strains == seqid.ls, "ep_distance"]) > 0) {
    epitope.dist[results.row, "epitope.dist.any.ls"] <- epitope.dist.any[epitope.dist.any$strains == seqid.ls, "ep_distance"]
  }

  # LATE COMERS:  the number of designated founders
  num.founders$num.founders.all[results.row] <- sum (lineages[, 2] == pubid.seq.tmp)
  num.founders$num.founders.tfl[results.row] <- length (grep ("TFL", lineages[lineages[, 2] == pubid.seq.tmp, 3]))

  # LATE COMERS:  subtype
  subtype[results.row] <- subtypes.table[subtypes.table$pubid == pubid.raw.tmp, "subtype.single"]

  # LATE COMERS:  country
  country[results.row] <- site.info[site.info$pub_id == pubid.raw.tmp, "country_iso"]
}

# it's best to handle our sitewise indicator variables separately
residues.bsites.indicators <- NULL

# TIER 2 (3):  binary indicators of residues found at CD4/VRC01 binding sites
for (site.num in 1:length (cd4.vrc01.bsites)) {

  # identify the site
  site.tmp <- cd4.vrc01.bsites[site.num]

  # get our observed residues
  residues.tmp <- sites.residues[[site.tmp]]
  residues.gaptext.tmp <- residues.tmp
  residues.gaptext.tmp[residues.gaptext.tmp == "-"] <- "gap"

  # sort out our alignment index
  index.study.tmp <- hxb2.map.study[hxb2.map.study[, 2] == site.tmp, 1]

  # create our temporary results object
  residues.bsites.indicators.tmp <- as.data.frame (matrix (NA, nrow=nrow (data), ncol=length (residues.tmp) * 3))
  names (residues.bsites.indicators.tmp) <- c (paste (paste (as.vector (sapply (paste ("hxb2", site.tmp, "is", sep="."), rep, length (residues.tmp))), residues.gaptext.tmp, sep="."), "mf", sep="."),
                                               paste (paste (as.vector (sapply (paste ("hxb2", site.tmp, "is", sep="."), rep, length (residues.tmp))), residues.gaptext.tmp, sep="."), "ms", sep="."),
                                               paste (paste (as.vector (sapply (paste ("hxb2", site.tmp, "is", sep="."), rep, length (residues.tmp))), residues.gaptext.tmp, sep="."), "ls", sep="."))
  
  # loop through each pubid
  for (pubid.num in 1:length (pubids)) {

    # set up some convenience variables
    pubid.raw.tmp <- pubids[pubid.num]
    pubid.seq.tmp <- paste0 ("V", gsub ("-", "_", pubid.raw.tmp))
    lineages.tmp <- lineages[lineages$pubid == pubid.seq.tmp, ]
    par.tmp <- par.scores[par.scores$seqname %in% lineages.tmp$seqid.aa, ]
    lineages.tmp <- merge (lineages.tmp, par.tmp[, c ("seqname", "pred.ic80")], by.x="seqid.aa", by.y="seqname")

    # determine the row of the results we want to modify
    results.row <- data$pub_id == pubid.raw.tmp

    # identify our representative sequences
    seqid.mf <- sample (lineages.tmp[lineages.tmp$seq.freq.aa == max (lineages.tmp$seq.freq.aa), "seqid.aa"], 1)
    seqid.ms <- sample (lineages.tmp[lineages.tmp$pred.ic80 == min (lineages.tmp$pred.ic80), "seqid.aa"], 1)
    seqid.ls <- sample (lineages.tmp[lineages.tmp$pred.ic80 == max (lineages.tmp$pred.ic80), "seqid.aa"], 1)

    # note the actual representative sequences
    seq.mf <- as.vector (seq.selected[[seq.id.lookup[[seqid.mf]]]])
    seq.ms <- as.vector (seq.selected[[seq.id.lookup[[seqid.ms]]]])
    seq.ls <- as.vector (seq.selected[[seq.id.lookup[[seqid.ls]]]])

    # capture our residues
    seq.residue.mf.tmp <- seq.mf[index.study.tmp]
    seq.residue.ms.tmp <- seq.ms[index.study.tmp]
    seq.residue.ls.tmp <- seq.ls[index.study.tmp]

    # determine and immortalize our indicator values
    seq.matches.tmp.mf <- seq.residue.mf.tmp == residues.tmp
    seq.matches.tmp.ms <- seq.residue.ms.tmp == residues.tmp
    seq.matches.tmp.ls <- seq.residue.ls.tmp == residues.tmp
    residues.bsites.indicators.tmp[results.row, ] <- as.numeric (c (seq.matches.tmp.mf, seq.matches.tmp.ms, seq.matches.tmp.ls))
  }

  # remove columns with zero counts
#  residues.bsites.indicators.tmp[, colSums (residues.bsites.indicators.tmp, na.rm=T) != 0]

  # build our dataset
  if (is.null (residues.bsites.indicators)) {
    residues.bsites.indicators <- residues.bsites.indicators.tmp
  } else {
    residues.bsites.indicators <- data.frame (residues.bsites.indicators, residues.bsites.indicators.tmp)
  }
}

# compile the RSA feature and add it to the data template
data.geo <- merge (data$pub_id, geo[, c ("pub_id", "country_iso")], by=1)
data.geo$is.rsa <- as.numeric (data.geo$country_iso == "ZA")
data <- data.frame (data[, 1:2], southAfrica=data.geo$is.rsa, data[, 3:ncol (data)])

# merge with the main dataset and export
#data.out <- cbind (data, residues.preselect, pngs, virgeom, residues.bsites, hdist)
data.out <- cbind (data, country, num.founders, subtype, residues.preselect, pngs, virgeom, residues.bsites.indicators, hdist, founders.sensitive.indicators, epitope.dist, glycosite.230)
outfilename <- paste0 ("amp_sieve_pooled_marks_final_v9.csv")
setwd ("~/Desktop")
write.csv (data.out, file=outfilename, row.names=F)

# break out sets by study and save
data.out.v703 <- data.out[data.out$protocol == "HVTN 703", ]
data.out.v704 <- data.out[data.out$protocol == "HVTN 704", ]
setwd (file.path (path.home, "out"))
write.csv (data.out.v703, file="amp_sieve_v703_marks_final_v9.csv", row.names=F)
write.csv (data.out.v704, file="amp_sieve_v704_marks_final_v9.csv", row.names=F)


# ---------------------------------------------------------------------------- #
#                                    - 30 -
# ---------------------------------------------------------------------------- #







