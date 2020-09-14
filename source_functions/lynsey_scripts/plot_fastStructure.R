setwd("/Volumes/mug01_n/taylorjerr/lwwvd/madtom/fastStructure/")

colors=c("#67001f","navy")

for(i in 1:2) {
	infile <- paste("MDTM.combined.all.filter.SNP.filtermaf.filtergeno.fs.", i, ".meanQ", sep="")
	outfile <- paste("MDTM.combined.all.filter.SNP.filtermaf.filtergeno.fs.", i, ".jpg", sep="")

	meanQ<-read.table(infile)
	breedNames<-read.table("fs_labels.txt", header=TRUE)

	jpeg(file=outfile, bg="white", width=720, height=460, units="px", pointsize=24)
	par(mai=c(2,2,0.5,0.5))
	barplot(t(as.matrix(meanQ)), width = 1, space = 0, beside = FALSE, col=colors, 
        xlab = "", ylab = "Ancestry", border = NA, xaxt="n", offset = 0)
	
	mtext(text=breedNames$breed, side=1, line=0, at=breedNames$n_for_label, las=3)
	abline(v=breedNames$running_total, lty=1, lwd=0.2)
	
	dev.off()
}

for(i in 1:3) {
  infile <- paste("MDTM.combined.all.filter.SNP.onlyneosho.filtermaf.filtergeno.fs.", i, ".meanQ", sep="")
  outfile <- paste("MDTM.combined.all.filter.SNP.onlyneosho.filtermaf.filtergeno.fs.", i, ".jpg", sep="")
  
  meanQ<-read.table(infile)
  breedNames<-read.table("fs_labels_onlyneosho.txt", header=TRUE)
  
  jpeg(file=outfile, bg="white", width=720, height=360, units="px", pointsize=14)
  par(mai=c(2,1,0.25,0.5))
  barplot(t(as.matrix(meanQ)), width = 1, space = 0, beside = FALSE, col=colors, 
          xlab = "", ylab = "Ancestry", border = NA, xaxt="n", offset = 0)
  
  mtext(text=breedNames$breed, side=1, line=0, at=breedNames$n_for_label, las=3)
  abline(v=breedNames$running_total, lty=1, lwd=0.2)
  
  dev.off()
}

for(i in 1:4) {
  infile <- paste("madtom_crtx_variants_passed_SNPs.fs.", i, ".meanQ", sep="")
  outfile <- paste("madtom_crtx_variants_passed_SNPs.fs.", i, ".jpg", sep="")
  
  meanQ<-read.table(infile)
  breedNames<-read.table("fs_labels.txt", header=TRUE)
  
  jpeg(file=outfile, bg="white", width=720, height=360, units="px", pointsize=14)
  par(mai=c(2,1,0.25,0.5))
  barplot(t(as.matrix(meanQ)), width = 1, space = 0, beside = FALSE, col=colors, 
          xlab = "", ylab = "Ancestry", border = NA, xaxt="n", offset = 0)
  
  mtext(text=breedNames$breed, side=1, line=0, at=breedNames$n_for_label, las=3)
  abline(v=breedNames$running_total, lty=1, lwd=0.2)
  
  dev.off()
}

for(i in 1:4) {
  infile <- paste("madtom_crtx_variants_passed_SNPs_onlyneosho_fs.", i, ".meanQ", sep="")
  outfile <- paste("madtom_crtx_variants_passed_SNPs_onlyneosho_fs.", i, ".jpg", sep="")
  
  meanQ<-read.table(infile)
  breedNames<-read.table("fs_labels_onlyneosho.txt", header=TRUE)
  
  jpeg(file=outfile, bg="white", width=720, height=360, units="px", pointsize=14)
  par(mai=c(2,1,0.25,0.5))
  barplot(t(as.matrix(meanQ)), width = 1, space = 0, beside = FALSE, col=colors, 
          xlab = "", ylab = "Ancestry", border = NA, xaxt="n", offset = 0)
  
  mtext(text=breedNames$breed, side=1, line=0, at=breedNames$n_for_label, las=3)
  abline(v=breedNames$running_total, lty=1, lwd=0.2)
  
  dev.off()
}