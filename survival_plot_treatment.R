treatment.survival.mCRC <- survfit(Surv(Diff, status) ~linje1kj,data=discovery.coxph)

ggsurvplot(treatment.survival.mCRC ,title="Survival by treatment categories", xlab="Time (Months)", ylab="Overall survival probability",
                font.main = 8, font.x =  8,font.y = 8,font.tickslab = 8,font.legend=8,pval.size = 2,pval.coord = c(10,0.3),
                size=0.4,legend = c(1,0.85),censor.size=2,break.time.by = 365.25*1, #tables.y.text = FALSE,
                pval= F, #"p=9.2e-06\nHR: 2.65\nCI: 1.43-4.9",
                palette = c("#1f78b4", "#9ecae1","#d7301f","#fc8d59","#bc80bd"),
                ggtheme = theme_classic(),risk.table = T, fontsize = 2 , tables.theme = theme_cleantable(),
                xscale="d_m",xlim=c(0,6*365), legend.title = "",
                legend.labs = c("5-FU",
                                "Oxaliplatin (Ox)",
                                "Irinotecan (Iri)",
                                "Ox+Iri",
                                "Immunotherapy")
)

p$table <- p$table +
  theme(plot.title = element_text(size = 8, color = "black", face = "bold"),axis.text.y = element_blank(),axis.ticks.y = element_line(colour = c("#bc80bd","#fc8d59", "#d7301f","#9ecae1","#1f78b4"),size=1),axis.ticks.length=unit(.4, "cm")) #+ theme(axis.ticks.y = element_blank())
