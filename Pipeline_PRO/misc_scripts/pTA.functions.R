adjacent.gene.coords <- function (fix.genes = NULL, bed.long = NULL, exon1 = NULL, bw.plus = NULL, 
    bw.minus = NULL, bp.bin = NULL, shift.up = NULL, dist.from.start = 50, 
    delta.tss = NULL, knot.div = 4, knot.thresh = 5, diff.tss = 1000, 
    pause.bp = 120, fname = "adjacentSplines.pdf") 
{
    TSS = c()
    pdf(fname)
    par(mfrow = c(3, 3))
    for (ii in 1:nrow(fix.genes)) {
        print(ii)
        first = fix.genes$upstream[ii]
        second = fix.genes$downstream[ii]
        ind.first = which(bed.long$gene == first)
        ind.second = which(bed.long$gene == second)
        bed.pair = bed.long[c(ind.first, ind.second), ]
        range = c(min(bed.pair$start), max(bed.pair$end))
        strand = bed.long$strand[ind.first]
        main = paste0(first, ", ", second, ", ", strand)
        if (strand == "+") {
            range[1] = range[1] - shift.up
            starts = seq(range[1], range[2] - bp.bin, by = bp.bin)
            ends = seq(range[1] + bp.bin, range[2], by = bp.bin)
            len = min(length(starts), length(ends))
            segments = data.frame(bed.pair$chr[1], starts[1:len], 
                ends[1:len], stringsAsFactors = FALSE)
            names(segments) = c("chr", "start", "end")
            segments$end[nrow(segments)] = range[2]
            counts = bed.region.bpQuery.bigWig(bw.plus, segments)
        }
        if (strand == "-") {
            range[2] = range[2] + shift.up
            ends = seq(range[2], range[1] + bp.bin, by = -bp.bin)
            starts = seq(range[2] - bp.bin, range[1], by = -bp.bin)
            len = min(length(starts), length(ends))
            segments = data.frame(bed.pair$chr[1], starts[1:len], 
                ends[1:len], stringsAsFactors = F)
            names(segments) = c("chr", "start", "end")
            segments$start[nrow(segments)] = range[1]
            counts = bed.region.bpQuery.bigWig(bw.minus, segments)
        }
        nknots = round(length(counts)/knot.div)
        if (nknots < knot.thresh) {
            nknots = knot.thresh
        }
        if (strand == "+") {
            spl = smooth.spline(x = segments$start, y = counts, 
                all.knots = FALSE, nknots = nknots)
            pks = findpeaks(spl$y, nups = 0)
            if(nrow(pks) > 1) {
                pks = pks[order(-pks[, 1]), ]
            }
            pk.coords = segments$start[pks[, 2]]
            ind.tss2 = which(abs(diff(pk.coords)) > diff.tss)[1] + 
                1
            tss1 = pk.coords[1]
            tss2 = pk.coords[ind.tss2]
            tss = c(tss1, tss2) - delta.tss
            spl = smooth.spline(x = segments$start, y = counts, 
                all.knots = FALSE, nknots = nknots)
            pred = predict(spl, segments$start)
            plot(segments$start, counts, main = main, xlab = "coords (bp)")
            points(pred$x, pred$y, type = "l", col = "red", xlab = "coords (bp)")
            abline(v = tss)
            new = cbind(c(first, second), c(min(tss), max(tss)))
        }
        if (strand == "-") {
            spl = smooth.spline(x = segments$end, y = counts, 
                all.knots = FALSE, nknots = nknots)
            pks = findpeaks(spl$y, nups = 0)
            if(nrow(pks) > 1) {   
                pks = pks[order(-pks[, 1]), ]
            }
            pk.coords = rev(segments$end)[pks[, 2]]
            ind.tss2 = which(abs(diff(pk.coords)) > diff.tss)[1] + 
                1
            tss1 = pk.coords[1]
            tss2 = pk.coords[ind.tss2]
            tss = c(tss1, tss2) + delta.tss
            spl = smooth.spline(x = segments$end, y = counts, 
                all.knots = FALSE, nknots = nknots)
            pred = predict(spl, segments$end)
            plot(segments$end, counts, main = main, xlab = "coords (bp)")
            points(pred$x, pred$y, type = "l", col = "red", xlab = "coords (bp)")
            abline(v = tss)
            new = cbind(c(second, first), c(min(tss), max(tss)))
        }
        new = as.data.frame(new, stringsAsFactors = FALSE)
        names(new) = c("gene", "TSS")
        TSS = rbind(TSS, new)
    }
    dev.off()
    TSS$TSS = as.numeric(TSS$TSS)
    new.bed = bed.long[bed.long$gene %in% c(fix.genes$upstream, 
        fix.genes$downstream), ]
    for (ii in 1:nrow(fix.genes)) {
        upgene = fix.genes$upstream[ii]
        dngene = fix.genes$downstream[ii]
        strand = bed.long$strand[bed.long$gene == upgene]
        pause.up = TSS$TSS[TSS$gene == upgene]
        pause.dn = TSS$TSS[TSS$gene == dngene]
        if (strand == "+") {
            tss.up = exon1$start[exon1$gene == upgene] %>% unique
            tss.dn = exon1$start[exon1$gene == dngene] %>% unique
            tss.up = tss.up[order(tss.up)]
            tss.dn = tss.dn[order(tss.dn)]
            dp_tss = pause.dn - tss.dn
            ind.dn = which(abs(dp_tss) == min(abs(dp_tss)))[1]
            tss.dn = tss.dn[ind.dn]
            tss.up = tss.up[tss.up < tss.dn]
            dp_tss = pause.up - tss.up
            ind.up = which(abs(dp_tss) == min(abs(dp_tss)))[1]
            tss.up = tss.up[ind.up]
            new.bed$start[new.bed$gene == upgene] = tss.up
            new.bed$start[new.bed$gene == dngene] = tss.dn
        }
        if (strand == "-") {
            tss.up = exon1$end[exon1$gene == upgene] %>% unique
            tss.dn = exon1$end[exon1$gene == dngene] %>% unique
            tss.up = tss.up[order(tss.up, decreasing = TRUE)]
            tss.dn = tss.dn[order(tss.dn, decreasing = TRUE)]
            dp_tss = tss.dn - pause.dn
            ind.dn = which(abs(dp_tss) == min(abs(dp_tss)))[1]
            tss.dn = tss.dn[ind.dn]
            tss.up = tss.up[tss.up > tss.dn]
            dp_tss = tss.up - pause.up
            ind.up = which(abs(dp_tss) == min(abs(dp_tss)))[1]
            tss.up = tss.up[ind.up]
            new.bed$end[new.bed$gene == upgene] = tss.up
            new.bed$end[new.bed$gene == dngene] = tss.dn
        }
    }
    for (ii in 1:nrow(fix.genes)) {
        upgene = fix.genes$upstream[ii]
        dngene = fix.genes$downstream[ii]
        strand = bed.long$strand[bed.long$gene == upgene]
        if (strand == "+") {
            tss.dn = new.bed$start[new.bed$gene == dngene]
            new.bed$end[new.bed$gene == upgene] = tss.dn - dist.from.start
        }
        if (strand == "-") {
            tss.dn = new.bed$end[new.bed$gene == dngene]
            new.bed$start[new.bed$gene == upgene] = tss.dn + 
                dist.from.start
        }
    }
    return(new.bed)
}
