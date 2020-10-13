
## Function 1: plotspectrum96 (COSMIC-like)

plotspectrum96 = function(m, pdfname="spectrum.pdf") {
    
    # Format context correctly
    if (is.null(m$context)) {
        stop("Trinucleotide context column ($context) missing from the input table")
    } else {
        m$context = paste(substr(m$context,1,1),substr(m$context,3,3),sep="-")
    }
    
    # Annotating the subs from the pyrimidine
    ntcomp = c(A="T",C="G",G="C",T="A")
    m$ref_py = m$ref
    m$mut_py = m$mut
    m$context_py = m$context
    for (j in 1:nrow(m)) {
        if (m$ref[j] %in% c("A","G")) {
            m$ref_py[j] = ntcomp[m$ref[j]]
            m$mut_py[j] = ntcomp[m$mut[j]]
            m$context_py[j] = paste(ntcomp[substr(m$context[j],3,3)],ntcomp[substr(m$context[j],1,1)], sep="-")
        }
    }
    
    # Calculating the frequencies
    freqs = table(paste(paste(m$ref_py,m$mut_py,sep=">"),m$context_py,sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
    
    dev.new(width=12,height=3)
    par(mar = c(4, 4.5, 2, 0.5))
    
    b = barplot(freqs_full, las=2, col=rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16), border=NA, space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", ylim=c(0,1.4*max(freqs_full)))
    rect(xleft=b[seq(1,81,by=16)]-0.5, xright=b[seq(16,96,by=16)]+0.5, ybottom=1.08*max(freqs_full), ytop=1.2*max(freqs_full), border=NA, col=c("dodgerblue","black","red","grey70","olivedrab3","plum2"))
    text(x=b[seq(8,96,by=16)]+1, y=1.3*max(freqs_full), labels=sub_vec)
    dev.copy(pdf, file=pdfname, width=12, height=3, useDingbats=F)
    #dev.off(); dev.off()
}

## Function 2: plotspectrum192
# plottype: (1) mirror plot, (2) two separate barplots

plotspectrum192 = function(m, pdfname="spectrum.pdf", trinfreqs=NULL, plottype=1) {
    
    # Format context correctly
    if (is.null(m$context)) {
        stop("Trinucleotide context column ($context) missing from the input table")
    } else {
        m$context = paste(substr(m$context,1,1),substr(m$context,3,3),sep="-")
    }
    
    # Annotating the subs from the pyrimidine
    ntcomp = c(A="T",C="G",G="C",T="A")
    m$ref_py = m$ref
    m$mut_py = m$mut
    m$context_py = m$context
    m$py = TRUE
    for (j in 1:nrow(m)) {
        if (m$ref[j] %in% c("A","G")) {
            m$ref_py[j] = ntcomp[m$ref[j]]
            m$mut_py[j] = ntcomp[m$mut[j]]
            m$context_py[j] = paste(ntcomp[substr(m$context[j],3,3)],ntcomp[substr(m$context[j],1,1)], sep="-")
            m$py[j] = FALSE
        }
    }
    
    # Calculating the frequencies PYR bases
    freqs = table(paste(paste(m$ref_py[m$py],m$mut_py[m$py],sep=">"),m$context_py[m$py],sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
    freqs_full_pyr = freqs_full
    
    # Calculating the frequencies PUR bases
    freqs = table(paste(paste(m$ref_py[!m$py],m$mut_py[!m$py],sep=">"),m$context_py[!m$py],sep=","))
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
    freqs_full_pur = freqs_full
    
    # Normalising for sequence composition
    if (!is.null(trinfreqs)) {
        tri = paste(substr(names(freqs_full_pyr),5,5), substr(names(freqs_full_pyr),1,1), substr(names(freqs_full_pyr),7,7), sep="")
        # Pyrimidine bases in the + strand
        freqs_full_pyr = freqs_full_pyr / trinfreqs[tri]
        # Purine bases in the + strand
        tricomp = sapply(tri, function(x) paste(ntcomp[rev(strsplit(x,split="")[[1]])],collapse=""))
        freqs_full_pur = freqs_full_pur / trinfreqs[tricomp]
        # Normalisation to sum to 1
        d = sum(freqs_full_pyr, freqs_full_pur)
        freqs_full_pyr = freqs_full_pyr/d
        freqs_full_pur = freqs_full_pur/d
        ylab = "Relative mutation rate"
    } else {
        ylab = "Number mutations"
    }
    
    if (plottype == 1) {
        
        dev.new(width=12,height=5)
        par(mar = c(4, 4.5, 2, 0.5))
        
        b = barplot(freqs_full_pyr, las=2, col=rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16), border=NA, space=1, cex.names=0.6, names.arg=xstr, ylab=ylab, ylim=c(-max(freqs_full_pyr,freqs_full_pur),1.3*max(freqs_full_pyr,freqs_full_pur)))
        barplot(-freqs_full_pur, add=T, col=rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16), border=NA, space=1, names.arg=NULL, xaxt='n', yaxt='n', ann=FALSE)
        rect(xleft=b[seq(1,81,by=16)]-0.5, xright=b[seq(16,96,by=16)]+0.5, ybottom=1.08*max(freqs_full_pyr), ytop=1.2*max(freqs_full_pyr), border=NA, col=c("dodgerblue","black","red","grey70","olivedrab3","plum2"))
        text(x=b[seq(8,96,by=16)]+1, y=1.3*max(freqs_full_pyr), labels=sub_vec)
        abline(h=0)
        dev.copy(pdf, file=pdfname, width=12, height=5, useDingbats=F)
        
    } else {
        
        dev.new(width=12,height=6)
        par(mfrow=c(2,1))
        par(mar = c(2.5, 4.5, 1.5, 0.5))
        
        maxy = max(freqs_full_pyr,freqs_full_pur)
        
        b = barplot(freqs_full_pyr, las=2, col=rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16), border=NA, space=1, cex.names=0.6, names.arg=xstr, ylab=ylab, ylim=c(0,1.4*maxy))
        rect(xleft=b[seq(1,81,by=16)]-0.5, xright=b[seq(16,96,by=16)]+0.5, ybottom=1.08*maxy, ytop=1.2*maxy, border=NA, col=c("dodgerblue","black","red","grey70","olivedrab3","plum2"))
        text(x=b[seq(8,96,by=16)]+1, y=1.3*maxy, labels=sub_vec)
        
        b = barplot(freqs_full_pur, las=2, col=rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16), border=NA, space=1, cex.names=0.6, names.arg=xstr, ylab=ylab, ylim=c(0,1.4*maxy))
        rect(xleft=b[seq(1,81,by=16)]-0.5, xright=b[seq(16,96,by=16)]+0.5, ybottom=1.08*maxy, ytop=1.2*maxy, border=NA, col=c("dodgerblue","black","red","grey70","olivedrab3","plum2"))
        text(x=b[seq(8,96,by=16)]+1, y=1.3*maxy, labels=sub_vec)
        dev.copy(pdf, file=pdfname, width=12, height=6, useDingbats=F)
        
    }
}

## Function 3: plotspectrum12

plotspectrum12 = function(m, pdfname="spectrum12.pdf", trinfreqs=NULL, conf.level=0.95, ymax=NULL) {
    
    nts = c("A","C","G","T")
    sub_types1 = c("C>A","C>G","C>T","T>A","T>C","T>G")
    sub_types2 = c("G>T","G>C","G>A","A>T","A>G","A>C")
    f = table(paste(m$ref[m$ref %in% nts & m$mut %in% nts], m$mut[m$ref %in% nts & m$mut %in% nts], sep=">"))
    f2 = rbind(f[sub_types1],f[sub_types2])
    f2[is.na(f2)] = 0
    
    # Confidence intervals 95%
    ci = sapply(c(f2), function(x) poisson.test(x, conf.level = conf.level)$conf.int)
    cilow = array(ci[1,], dim=c(2,6))
    cihigh = array(ci[2,], dim=c(2,6))
    
    # Sequence composition normalisation
    if (!is.null(trinfreqs)) {
        fseq = table(rep(substr(names(trinfreq),2,2), times=trinfreq))[nts]
        fseq[is.na(fseq)] = 0
        names(fseq) = nts
        fseq = array(fseq[c(rep(c("C","G"),3),rep(c("T","A"),3))], dim=c(2,6))
        
        cilow = cilow / fseq; cilow = cilow / sum(f2/fseq)
        cihigh = cihigh / fseq; cihigh = cihigh / sum(f2/fseq)
        f2 = f2 / fseq; f2 = f2 / sum(f2)
    }
    
    # Plot
    dev.new(width=4,height=5)
    if (is.null(ymax)) {
        h = barplot(f2, beside = T, border=NA, las=2, col=c("grey20","grey70"), names.arg=rep("",6), ylim=c(0,max(cihigh)), 
                    xlab="Mutation types (plus strand)", ylab=c("Mutations observed","Relative mutation rate")[(!is.null(trinfreqs))+1])
    } else {
        h = barplot(f2, beside = T, border=NA, las=2, col=c("grey20","grey70"), names.arg=rep("",6), ylim=c(0,ymax), 
                    xlab="Mutation types (plus strand)", ylab=c("Mutations observed","Relative mutation rate")[(!is.null(trinfreqs))+1])
    }
     axis(1, at=h[1,], tick=F, labels=sub_types1, las=2)
    axis(1, at=h[2,], tick=F, labels=sub_types2, las=2)
    segments(x0=h, y0=cilow, y1=cihigh)
    dev.copy(pdf, file=pdfname, width=4, height=5, useDingbats=F)
    
    # Output the substitution table
    return(f2)
}

