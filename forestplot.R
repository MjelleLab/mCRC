mCRC.ss.extended.discovery.order.ss.order #metadata file
mCRC.matrix.survival.inclusion.dge.cpm.discovery.order.t #transposed count file, cpm normalized

discovery.coxph <- as.data.frame(merge(mCRC.ss.extended.discovery.order.ss.order,mCRC.matrix.survival.inclusion.dge.cpm.discovery.order.t,by.y="row.names",by.x="ID2"))
discovery.coxph <- discovery.coxph[,grep("hsa|Diff|^status|kjonn|age|cea_base|crp_base|linje1kj|msi_mmr",colnames(discovery.coxph))]
discovery.coxph$status <- gsub("2","0",discovery.coxph$status)
discovery.coxph$status  <- as.numeric(discovery.coxph$status)

coxph.Inclusion <- capture.output(for(i  in colnames(discovery.coxph)){
  print(summary(coxph(as.formula(paste0("Surv(Diff, status)~", i,"+kjonn","+age","+linje1kj","+msi_mmr")),  data=as.data.frame(discovery.coxph))))
})



out_coxph.discovery_sRNA_Biobank1 <- as.matrix(as.list(coxph.Inclusion)[grep("hsa", as.list(coxph.Inclusion))])  
out_coxph.discovery_sRNA_Biobank1 <- as.data.frame(gsub("\\s+", " ", out_coxph.discovery_sRNA_Biobank1)  )      #Convert spaces to columns
toDelete <- seq(1, nrow(out_coxph.discovery_sRNA_Biobank1), 2)
toDelete2 <- seq(0, nrow(out_coxph.discovery_sRNA_Biobank1), 2)

out_coxph.discovery_sRNA_Biobank1.1 <- as.data.frame(out_coxph.discovery_sRNA_Biobank1[ toDelete ,])
out_coxph.discovery_sRNA_Biobank1.2 <- as.data.frame(out_coxph.discovery_sRNA_Biobank1[ toDelete2 ,])

out_coxph.discovery_sRNA_Biobank1.1 <- splitstackshape::cSplit(out_coxph.discovery_sRNA_Biobank1.1, 'out_coxph.discovery_sRNA_Biobank1[toDelete, ]', sep = ' ')
out_coxph.discovery_sRNA_Biobank1.2 <- splitstackshape::cSplit(out_coxph.discovery_sRNA_Biobank1.2, 'out_coxph.discovery_sRNA_Biobank1[toDelete2, ]', sep = ' ')

colnames(out_coxph.discovery_sRNA_Biobank1.1)  <- c("Parameter","coef","hz","se_coef","z","p_val","Sign")
colnames(out_coxph.discovery_sRNA_Biobank1.2)  <- c("Parameter","coef","-hz","lower","upper")
out_coxph.discovery_sRNA_Biobank1.out <- as.data.frame(cbind(out_coxph.discovery_sRNA_Biobank1.1,out_coxph.discovery_sRNA_Biobank1.2))
out_coxph.discovery_sRNA_Biobank1.out <- out_coxph.discovery_sRNA_Biobank1.out[order(out_coxph.discovery_sRNA_Biobank1.out$p_val, decreasing = F),]
out_coxph.discovery_sRNA_Biobank1.out
out_coxph.discovery_sRNA_Biobank1.out$BH_Adj_Pval <- p.adjust(out_coxph.discovery_sRNA_Biobank1.out$p_val,method = "BH")
out_coxph.discovery_sRNA_Biobank1.out.sign <- out_coxph.discovery_sRNA_Biobank1.out[out_coxph.discovery_sRNA_Biobank1.out$BH_Adj_Pval<0.05,]
out_coxph.discovery_sRNA_Biobank1.out.sign <- out_coxph.discovery_sRNA_Biobank1.out.sign[order(out_coxph.discovery_sRNA_Biobank1.out.sign$BH_Adj_Pval,decreasing = F),]
head(out_coxph.discovery_sRNA_Biobank1.out.sign,10)
out_coxph.discovery_sRNA_Biobank1.out.sign$Sign_Adj <- out_coxph.discovery_sRNA_Biobank1.out.sign$BH_Adj_Pval<0.05
head(out_coxph.discovery_sRNA_Biobank1.out.sign,30)

coxph.in.consistent.in.discovery <- out_coxph.discovery_sRNA_Biobank1.out.sign

Covariate.mirna <- gsub("hsa-","",gsub("_","-",coxph.in.consistent.in.discovery$Parameter))
HR.disc.miRNA <- coxph.in.consistent.in.discovery$hz
lower.disc.miRNA <- coxph.in.consistent.in.discovery$lower
upper.disc.miRNA <- coxph.in.consistent.in.discovery$upper
pvalue.disc.miRNA <- ifelse(coxph.in.consistent.in.discovery$BH_Adj_Pval <0.005, format(as.numeric(coxph.in.consistent.in.discovery$BH_Adj_Pval),scientific = T,digits = 3), round(coxph.in.consistent.in.discovery$BH_Adj_Pval,4))





dummydata.miRNA.inclusion <- data.frame(Covariate.mirna, 
                                        HR.disc.miRNA, 
                                        lower.disc.miRNA, 
                                        upper.disc.miRNA, 
                                        pvalue.disc.miRNA
                                        )




tabletextdummy3.miRNA.inclusion <- cbind(c("MicroRNA",NA,Covariate.mirna),
                                         cbind(c("P-value\t\t\t\t\tHR\t\t\t\t\t\t\t\t\t\t\t\tCI",NA,paste(paste("P:",format(as.numeric(pvalue.disc.miRNA),scientific=T),sep=""),
                                                                      paste(as.character(paste("\t\tHR:",sprintf("%.2f",round(HR.disc.miRNA,2))),sep=""),paste(paste("\t\t\t[",paste(round(lower.disc.miRNA,2),round(upper.disc.miRNA,2),sep=","),sep=""),"]",sep=""),sep=" "),
                                                                      sep=" "))))




forestplot(mean= c(NA,NA, dummydata.miRNA.inclusion$HR.disc.miRNA),
           lower = c(NA,NA,dummydata.miRNA.inclusion$lower.disc.miRNA), 
           upper = c(NA,NA,dummydata.miRNA.inclusion$upper.disc.miRNA),
           tabletextdummy3.miRNA.inclusion, 
           new_page = TRUE,
           align="l",
           clip = c(0.1,5), 
           #lineheight = unit(10,"mm"),
           line.margin = .1,
           xlog = TRUE, xlab = "HR with 95% CI", 
           col = fpColors(box = c("#d73027", "#4575b4"), 
                          lines = c("grey80", "grey80")),
           fn.ci_norm = "fpDrawNormalCI", #           fn.ci_norm = c(fpDrawNormalCI, fpDrawDiamondCI)
           is.summary = c(TRUE,rep(FALSE,28)), 
           graph.pos = 2,
           boxsize = 0.4, 
           xticks = c(0.5, 1, round(max(c(dummydata.miRNA.inclusion$upper.disc.miRNA,na.rm = T)))),
           vertices = TRUE,
           txt_gp=fpTxtGp(xlab  = gpar(cex = 0.75),
                          cex=0.65,
                          summary=gpar(cex = 0.75),
                          ticks = gpar(fontfamily = "", cex = 0.55),
                          label = gpar(fontfamily = "Helvetica"))
           )
