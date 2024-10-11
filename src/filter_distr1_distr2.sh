#!/bin/bash
# author : sjn

# sample distribution has relprobs of [0.1, 0.2, 0.1, 0.05, 0.05, 0.3, 0.05, 0.05, 0.05, 0.05]
# real distribution has counts of [20k, 4k, 4k, 900, 900, 900, 900, 900, 4k, 4k]
# assign 0.3 (highest relprob) from sd to 900 count from corresponding rd bin
# all other relprobs (ex: 0.05) can be used as in: 0.3/0.05 = 900/x => solve for x (target count from rd)
#
# shuffle all real data tags.
# for each interval/break-bin, keep count as you traverse file
# stop output for a particular interval/break-bin when you reach its target count value

# on this initial code write:
# if a particular interval/break-bin cannot meet target count value, so be it - output everything for that interval/break-bin
# more generally, reduce the count (900 in example above) so that every interval can meet its recalculated target count

# the subsampling from samtools does not re-order anything, and so my assumptions were incorrect - some things below are
#  not doing what I thought and/or aren't needed.  For example, sampling with fraction=1 just returns a copy of the original

if [[ $# != 3 ]]; then
  printf "Expect $0 <sample-distr.ft.bam> <input.ft.bam> <output-filtered.ft.bam>\n"
  exit 1
fi

set -exuo pipefail

sample_distribution=$1 # .ft.bam file
inp=$2 # .ft.bam file
out=$3

if [[ ! -s ${inp} ]]; then
  printf "Problem finding 1 file: %s\n" ${inp}
  exit 1
fi

TMPDIR="/tmp"

tmpd=${TMPDIR}/$(whoami)/$$
rm -rf ${tmpd}
mkdir -p ${tmpd}
mkdir -p $(dirname "${out}")

ftype=m6a

ft extract ${sample_distribution} --all - \
  | cutnm total_m6a_bp,total_AT_bp \
 >${tmpd}/sample_distr.${ftype}

cat ${inp} \
  | ft extract --all - \
  | cutnm total_m6a_bp,total_AT_bp \
  >${tmpd}/input.${ftype}

R --no-save --quiet <<__R__
  minx <- 0.00
  maxx <- 1
  sample_distr <- read.table("$tmpd/sample_distr.$ftype", header=TRUE, sep="\t", row.names=NULL)
  real_data <- read.table("$tmpd/input.$ftype", header=TRUE, sep="\t", row.names=NULL)
  ss <- sample_distr[1]/sample_distr[2]
  dd <- real_data[1]/real_data[2]
  colnames(ss) <- NULL
  colnames(dd) <- NULL
  sample_distr <- unlist(ss)
  real_data <- unlist(dd)

  delta <- 0.001
  brks <- seq(0, maxx+delta, delta)
  h <- hist(sample_distr, breaks=brks, plot=FALSE)
  j <- hist(real_data, breaks=brks, plot=FALSE)

  rel_probs <- h[["counts"]]/length(sample_distr)
  # adding offset from 5 to keep anomalous things near 0/1 from hurting things
  mx_index_rel_probs <- which.max(rel_probs[5:length(rel_probs)]) # match to actual count in real_data and calculate all counts from this
  root_prob <- rel_probs[mx_index_rel_probs]
  root_count <- j[["counts"]][mx_index_rel_probs]
  fi <- findInterval(real_data, brks) # gives index of brks for each real_data

  rd <- cbind(real_data, fi)
  colnames(rd) <- c("real_data", "rd_interval")
  rd <- as.data.frame(rd)

  targets <- floor(root_count * rel_probs / root_prob)

  ## This SLOW loop shows what's going on below with ave() that cheatgpt showed as an equiv
  # vals <- c()
  # for (idx in 1:(dim(rd)[1])) {
  #   targets[rd[idx,2]] <- targets[rd[idx,2]]-1
  #   vals <- c(vals, targets[rd[idx,2]]>=0)
  # }

  # Use ave() to apply the decrement per group, replicating the behavior of a loop
  # Create a helper function to decrement each index correctly
  target_indices <- rd[, 2]
  targets[target_indices] <- ave(target_indices, target_indices, FUN=function(x) targets[x[1]] - seq_along(x))
  vals <- targets[target_indices] >= 0

  target_indices <- rd[, 2]
  targets[target_indices] <- targets[target_indices] - 1
  vals <- targets[target_indices] >= 0

  vals <- as.data.frame(vals*1)
  write.table(vals, file="$tmpd/output.$ftype", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  f <- cbind(real_data, vals)
  colnames(f) <- c("d", "v")
  f <- unname(unlist(subset(f, v==1, select=d)))
  f <- as.numeric(f)
  pdf(paste("$out", ".pdf", sep=""))
  par(mfrow=c(3,1))
  hist(sample_distr, breaks=brks, xlim=c(0,0.3))
  hist(real_data, breaks=brks, xlim=c(0,0.3))
  hist(f, breaks=brks, xlim=c(0,0.3), main="real_data under sample_distr")
  dev.off()
__R__

samtools view -H ${inp} \
 >${tmpd}/.header

samtools view ${inp} \
  | paste ${tmpd}/output.$ftype - \
  | awk '$1 > 0' \
  | cut -f2- \
  | cat ${tmpd}/.header - \
  | samtools view -h -b \
  >${out}

pbindex ${out}

rm -rf ${tmpd}

exit 0
