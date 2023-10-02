library("RColorBrewer")
library("ggplot2")
library(tidyverse)
library(dplyr)
library(ggalt)
library(ggpubr)
library(scales)
library(pheatmap)
library(ComplexHeatmap)
library(ggformula)
library(gridExtra)

cal_z_score <- function(x){
    (x-mean(x)) / sd(x)
}
logfnxn = function(x){
    if(x==0){
        x = 0.1
    }
    x = as.numeric(as.character(x))
    log10(x)
}
pal_sampletype <- c("Blood"="darkorange1","BAL"="royalblue2")

pal_deadoralive <- c("Alive"="#776bcd","Dead"="firebrick1")

pal_protocol <-  c("Control"="#22BB3B","COVID"="#C40F5B")

####################################################################################
#Set up data objects for line plots and box plots
####################################################################################
useForLinesOrBoxplots = function(only_baseline, include_controls){
    # File: sample_ig_levels.csv [rows contain antibody isotype + sample type + sample ID + Ab levels for each subfraction/antibody analyzed]
    dat_IGs <- read.csv("sample_ig_levels.csv", sep=",")
    dat_IGs <- dat_IGs %>% pivot_longer(    #Re-format file from wide to long
        cols=-c(dilution_to_use,protocol,AB,plate,Sample_type,yl_sample_id,SamplePlasma,Dilution),
        names_to="AB_name",
        values_to="AB_value"
        ) %>% filter(dilution_to_use==1)
    # File: sample_mappingfile.csv.csv [rows contain sample ID and then sample/subject metadata]
    dat_smpls <- read.csv("COVID_sample_mappingfile.csv", sep=",")

    #Set include_controls=TRUE if analysis compares cases to controls
    if(include_controls){   
        copd_mapping <- read.csv("COPD_sample_mappingfile.csv", sep=",")        
        dat_smpls$protocol = 'COVID'
        copd_mapping$protocol="Control"
        dat_smpls <- rbind(copd_mapping, dat_smpls)
    }

    #Set only_baseline=TRUE if analysis should only include the baseline (aka first) sample for given subject
    if(only_baseline){
        dat_smpls <- dat_smpls %>% group_by(sample_type, study_id) %>% slice_min(adm_to_sample_days) %>% ungroup()
    }

    rownames(dat_smpls) <- dat_smpls$yl_sample_id

    #Create dat object by merging the object with sample metadata to the one containing the antibody levels
    if(include_controls){
        dat <- merge(dat_smpls, dat_IGs%>%select(-c("protocol")), by=c("yl_sample_id"))
    }else{
        dat <- merge(dat_smpls, dat_IGs, by=c("yl_sample_id"))    
    }
    dat$AB <- factor(dat$AB,levels=c("IgA","IgM",'IgG'))

    #AB level should never be zero, so add 0.1 to any that are (0.1 is the lowest detected level)
    dat <- dat %>% mutate(AB_value_plus1=case_when(
        AB_value == 0 ~ 0.1,
        TRUE ~ AB_value

    ))
    return(dat)
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Create figure: PAIRED BAL/BLOOD BOXPLOTS
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#dat contains a row for every unique sample/IG isotype/AB epitope combination
dat <- useForLinesOrBoxplots(TRUE,FALSE) #TRUE to only include baseline samples, FALSE to not include controls

#Find only the samples that are paired (BAL and blood collected from same subject on same day)
dat_paired_wide <- dat %>% pivot_wider(
    id_cols=c("study_id","adm_to_sample_days","AB","AB_name"),
    names_from=c(Sample_type),
    values_from=c(yl_sample_id),
    names_prefix="yl_sample_id_"
    ) %>% filter(!is.na(yl_sample_id_BAL), !is.na(yl_sample_id_Plasma))
# Don't need this because already only have baselines: dat_paired_wide <- dat_paired_wide %>% group_by(study_id,AB,AB_name) %>% slice_min(adm_to_sample_days) %>% ungroup()

#dat_paired_post contains a row for every sample that is part of a pair described above
dat_paired_post <- dat %>% filter((yl_sample_id %in% dat_paired_wide$yl_sample_id_BAL) | (yl_sample_id %in% dat_paired_wide$yl_sample_id_Plasma))

dat_paired_post$sample_type = factor(dat_paired_post$sample_type,levels=c("Blood","BAL"))
dat_paired_post$AB_name = factor(dat_paired_post$AB_name,levels=c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid","TT","LukS"))
txt_size=22

a<-ggplot(dat_paired_post %>% filter(AB_name=='Spike'), aes(x = sample_type, y = AB_value_plus1)) +
    geom_boxplot(aes(fill = sample_type,color=sample_type), alpha = 0.3,width=0.7)  +geom_line(aes(group = study_id),color="gray60")+geom_point(aes(col = sample_type))+
    scale_color_manual(values=pal_sampletype,name="")+scale_fill_manual(values=pal_sampletype,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("Spike")+labs(y='Log Immunoglobulin level')+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size,margin=margin(r=10)),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=2),panel.spacing=unit(0.2,"lines"))+#t,r,b,l
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10,paired=TRUE)+
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7,paired=TRUE)

b<-ggplot(dat_paired_post %>% filter(AB_name=='SpikeNTD'), aes(x = sample_type, y = AB_value_plus1)) +
    geom_boxplot(aes(fill = sample_type,color=sample_type), alpha = 0.3,width=0.7)  +geom_line(aes(group = study_id),color="gray60")+geom_point(aes(col = sample_type))+
    scale_color_manual(values=pal_sampletype,name="")+scale_fill_manual(values=pal_sampletype,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("SpikeNTD")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=2),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10,paired=TRUE)+
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7,paired=TRUE)

c <- ggplot(dat_paired_post %>% filter(AB_name=='SpikeRBD'), aes(x = sample_type, y = AB_value_plus1)) +
    geom_boxplot(aes(fill = sample_type,color=sample_type), alpha = 0.3,width=0.7) +geom_line(aes(group = study_id),color="gray60")+geom_point(aes(col = sample_type)) +
    scale_color_manual(values=pal_sampletype,name="")+scale_fill_manual(values=pal_sampletype,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("SpikeRBD")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=2),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10,paired=TRUE)+
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7,paired=TRUE)

d <- ggplot(dat_paired_post %>% filter(AB_name=='Nucleocapsid'), aes(x = sample_type, y = AB_value_plus1)) +
    geom_boxplot(aes(fill = sample_type,color=sample_type), alpha = 0.3,width=0.7) +geom_line(aes(group = study_id),color="gray60")+geom_point(aes(col = sample_type)) +
    scale_color_manual(values=pal_sampletype,name="")+scale_fill_manual(values=pal_sampletype,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("Nucleocapsid")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=2),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10,paired=TRUE)+
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7,paired=TRUE)

a <- a + scale_y_log10(limits=c(5,100000), labels=trans_format('log10', math_format(10^.x)))
b <- b + scale_y_log10(limits=c(5,100000), labels=trans_format('log10', math_format(10^.x)))
c <- c + scale_y_log10(limits=c(5,100000), labels=trans_format('log10', math_format(10^.x)))
d <- d + scale_y_log10(limits=c(5,100000), labels=trans_format('log10', math_format(10^.x)))

pdf("pairedboxplot_topo_bl_covidrelated_wpairedstats.pdf",width=12.45, height=4.4)
# pdf("pairedboxplot_topo_bl_covidrelated_notpairedstats.pdf",width=12.45, height=4.4)
gridExtra::grid.arrange(
    a,b,c,d,
    ncol=4,
    widths=c(3.45,3,3,3),
    heights=c(4.4))
dev.off()


a<-ggplot(dat_paired_post %>% filter(AB_name=='TT'), aes(x = sample_type, y = AB_value_plus1)) +
    geom_boxplot(aes(fill = sample_type,color=sample_type), alpha = 0.3,width=0.7)+geom_line(aes(group = study_id),color="gray60") +geom_point(aes(col = sample_type)) +
    scale_color_manual(values=pal_sampletype,name="")+scale_fill_manual(values=pal_sampletype,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("TT")+labs(y='Log Immunoglobulin level')+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size,margin=margin(r=10)),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=2),panel.spacing=unit(0.2,"lines"))+#t,r,b,l
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10)+
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7)

b<-ggplot(dat_paired_post %>% filter(AB_name=='LukS'), aes(x = sample_type, y = AB_value_plus1)) +
    geom_boxplot(aes(fill = sample_type,color=sample_type), alpha = 0.3,width=0.7)  +geom_line(aes(group = study_id),color="gray60")+geom_point(aes(col = sample_type))+
    scale_color_manual(values=pal_sampletype,name="")+scale_fill_manual(values=pal_sampletype,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("LukS")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=2),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10)+#,paired=TRUE
    stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7)

a <- a + scale_y_log10(limits=c(2,30000), labels=trans_format('log10', math_format(10^.x)))
b <- b + scale_y_log10(limits=c(2,30000), labels=trans_format('log10', math_format(10^.x)))

# pdf("pairedboxplot_topo_bl_notcovidrelated_wpairedstats.pdf",width=6.45, height=4.4)
pdf("pairedboxplot_topo_bl_notcovidrelated_notpairedstats.pdf",width=6.45, height=4.4)
gridExtra::grid.arrange(
    a,b,
    ncol=2,
    widths=c(3.45,3),
    heights=c(4.4))
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Create figure: TOPOGRAPHIC HEATMAP
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


OG_file <- read.csv("sample_ig_levels.csv", sep=,)%>% filter(dilution_to_use==1)
OG_file <- OG_file %>% filter(protocol=='COVID')
OG_file <- OG_file %>% select("AB","yl_sample_id","Spike","SpikeNTD","SpikeRBD","Nucleocapsid","S435.462","S435.462.N439K.Y543F","S524.560","S564.584","S564.584.UK","S605.636.ADE","S802.827","S888.909","S1184.1209","LukS","TT")

OG_file <- OG_file %>% pivot_longer(    #Re-format file from wide to long
        cols=c(Spike,SpikeNTD,SpikeRBD,Nucleocapsid,S435.462,S435.462.N439K.Y543F,S524.560,S564.584,S564.584.UK,S605.636.ADE,S802.827,S888.909,S1184.1209,LukS,TT),
        names_to="AB_name",
        values_to="AB_value"
        ) 
OG_file$AB_name = paste(OG_file$AB_name,"_",OG_file$AB,sep="")
OG_file <- as.data.frame(OG_file %>% select(-c("AB")) %>% pivot_wider(id_cols=yl_sample_id,names_from=AB_name,values_from=AB_value))
rownames(OG_file) <- OG_file$yl_sample_id
dat_IGs <- OG_file%>%select(-c("yl_sample_id"))

# File: sample_mappingfile.csv.csv [rows contain sample ID and then sample/subject metadata]
dat_smpls <- read.csv("COVID_sample_mappingfile.csv", sep=",")

#****************************#
#*** OPTION OPTION OPTION ***#
#Use this if only want to use baseline
dat_smpls <- dat_smpls %>% group_by(sample_type, study_id) %>% slice_min(adm_to_sample_days) %>% ungroup()
#****************************#

needed <- which(rownames(dat_IGs) %in% dat_smpls$yl_sample_id)
dat_IGs <- dat_IGs[needed,]

#****************************#
#*** OPTION OPTION OPTION ***#
#Choose scoring system

#Get Z Score
# data_norm <- t(apply(dat_IGs,1,cal_z_score))
# dat_IGs <- data_norm

#Get Log10 - I add 1 to everything to deal with the values of 0
dat_IGs <- apply(dat_IGs, c(1,2), logfnxn)

#Get scaled 
#  dat_IGs_noNA_scale = scale(dat_IGs)
# dat_IGs <- dat_IGs_noNA_scale
 #****************************#

#Make sure ordered okay
dat_smpls <- dat_smpls[which(dat_smpls$yl_sample_id %in% rownames(dat_IGs)),]
dat_smpls <- dat_smpls[order(dat_smpls$yl_sample_id),]
dat_IGs <-dat_IGs[order(rownames(dat_IGs)),]

#Transform so rows are "Spike_IgA, etc" and cols are samples
dat_IGs_t <- as.data.frame(t(as.matrix(dat_IGs))) #All chosen IG

annon_colors= list(sample_type=c(BAL="royalblue2",Blood="darkorange1")) 
   
annotation_col_sampletype <- data.frame(
    sample_type = dat_smpls$sample_type,
    row.names = colnames(dat_IGs_t)
)

sampletype_cols_list <-  ifelse(
    grepl('BAL',colnames(dat_IGs_t))==TRUE,"BAL",
    ifelse(grepl('Plasma',rownames(dat_IGs_t)),"Blood","Blood")
)


pdf("heatmap_topo_bl_log10_nolegend_nosubjectID.pdf",height = 16, width = 20) 
    ComplexHeatmap::pheatmap(dat_IGs_t,annotation_col=annotation_col_sampletype,annotation_colors=annon_colors,cluster_rows=FALSE,gaps_row=c(15,30),legend=FALSE,annotation_legend=FALSE,show_colnames=FALSE,
    column_title = "Baseline samples", column_title_side = "bottom",
    fontsize=20,
    column_title_gp = gpar(fontsize = 25),
    treeheight_col=unit(2.6,"cm")
    )
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Create figure: DEAD VS ALIVE 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#dat contains a row for every unique sample/IG isotype/AB epitope combination
dat <- useForLinesOrBoxplots(TRUE,FALSE) #TRUE to only include baseline samples, FALSE to not include controls

dat$AB_name = factor(dat$AB_name,levels=c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid","TT","LukS"))
dat$dead_or_alive = factor(dat$dead_or_alive,levels=c("Alive","Dead"))
dat <- dat %>% mutate(subfraction_category=case_when(
    AB_name == 'Spike' | AB_name == 'SpikeNTD' | AB_name == 'SpikeRBD' | AB_name == 'Nucleocapsid' ~ 'covid_related',
    AB_name == 'S435-462' | AB_name == 'S435-462#N439K#Y543F' | AB_name == 'S524-560' | AB_name == 'S564-584' | AB_name == 'S564-584 UK' | AB_name == 'S605-636#ADE' | AB_name == 'S802-827' | AB_name == 'S888-909' | AB_name == 'S1184-1209' ~ 'spike_subfxn',
    AB_name != 'Nucleosome' & AB_name != 'SpikeS2' ~ 'non_covid_related', 
    TRUE ~ 'misc_unused'
))

## ## ## BAL ONES ## ## ##
    dat <- dat %>% filter(sample_type=='BAL')
## ## ## BLOOD ONES ## ## ##
    dat <- dat %>% filter(sample_type=='Blood')

    ### ### Create the baseline boxplots ### ###
    txt_size=22
    a<- ggplot(dat %>% filter(AB=='IgA') %>% filter(AB_name != 'NA'), aes(x=AB_name, y=AB_value_plus1,color=factor(dead_or_alive)))+
       geom_boxplot(aes(fill=dead_or_alive,color=dead_or_alive), width=0.5,position=position_dodge(0.7),outlier.shape = NA,alpha=0.3)+
       geom_point(position=position_jitterdodge(dodge.width=0.5),size=0.8,alpha=0.8)+
        scale_color_manual(values=pal_deadoralive,name="")+scale_fill_manual(values=pal_deadoralive,name="")+
        facet_grid(cols=vars(subfraction_category),scales="free_x",space="free")+
        ggtitle("IgA")+labs(y='Log Immunoglobulin level')+
        theme_pubr()+theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,1,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_text(angle = 45, hjust=1,size=txt_size,color="black"),axis.text.y=element_text(size=txt_size,,color="black"),axis.title.y=element_text(size=txt_size,margin=margin(r=10),vjust=1.5),axis.ticks.x=element_blank(),strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_blank(),panel.spacing=unit(0.5,"lines"))+
        stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10)+
        stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7)
    #BAL
        # a <- a + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))
    #Blood
        a <- a + scale_y_log10(limits=c(20,100000), labels=trans_format('log10', math_format(10^.x)))

    b<- ggplot(dat %>% filter(AB=='IgM') %>% filter(AB_name != 'NA'), aes(x=AB_name, y=AB_value_plus1,color=factor(dead_or_alive)))+
       geom_boxplot(aes(fill=dead_or_alive,color=dead_or_alive), width=0.5,position=position_dodge(0.7),outlier.shape = NA,alpha=0.3)+
       geom_point(position=position_jitterdodge(dodge.width=0.5),size=0.8,alpha=0.8)+
        scale_color_manual(values=pal_deadoralive,name="")+scale_fill_manual(values=pal_deadoralive,name="")+
        facet_grid(cols=vars(subfraction_category),scales="free_x",space="free")+
        ggtitle("IgM")+labs(y='')+
        theme_pubr()+theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,1,5.5,5.5,"pt"))+theme(axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_text(angle = 45, hjust=1,size=txt_size,color="black"),axis.text.y=element_text(size=txt_size,,color="black"),axis.title.y=element_text(size=txt_size,margin=margin(r=10),vjust=1.5),axis.ticks.x=element_blank(),strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_blank(),panel.spacing=unit(0.5,"lines"))+
        stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10,paired=TRUE)+
        stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7,paired=TRUE)
    #BAL
        # b <- b + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))
    # Blood
        b <- b + scale_y_log10(limits=c(20,100000), labels=trans_format('log10', math_format(10^.x)))

    c<- ggplot(dat %>% filter(AB=='IgG') %>% filter(AB_name != 'NA'), aes(x=AB_name, y=AB_value_plus1,color=factor(dead_or_alive)))+
       geom_boxplot(aes(fill=dead_or_alive,color=dead_or_alive), width=0.5,position=position_dodge(0.7),outlier.shape = NA,alpha=0.3)+
       geom_point(position=position_jitterdodge(dodge.width=0.5),size=0.8,alpha=0.8)+
        scale_color_manual(values=pal_deadoralive,name="")+scale_fill_manual(values=pal_deadoralive,name="")+
        facet_grid(cols=vars(subfraction_category),scales="free_x",space="free")+
        ggtitle("IgG")+labs(y='')+
        theme_pubr()+theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"))+theme(axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_text(angle = 45, hjust=1,size=txt_size,color="black"),axis.text.y=element_text(size=txt_size,,color="black"),axis.title.y=element_text(size=txt_size,margin=margin(r=10),vjust=1.5),axis.ticks.x=element_blank(),strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_blank(),panel.spacing=unit(0.5,"lines"))+
        stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10,paired=TRUE)+
        stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.04,size=7,paired=TRUE)
    #BAL
        # c <- c + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))
    #Blood
        c <- c + scale_y_log10(limits=c(20,100000), labels=trans_format('log10', math_format(10^.x)))

    # ggsave(filename=paste0("boxplot_outcome_bl_bal.pdf"),
    ggsave(filename=paste0("boxplot_outcome_bl_blood.pdf"),
       ggarrange(a,b,c,ncol=3,nrow=1),
    height=6,width=18)


### ### Create the line plots for BAL ### ###

dat <- useForLinesOrBoxoplots(FALSE,FALSE)
dat$AB_name = factor(dat$AB_name,levels=c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid"))
dat$dead_or_alive = factor(dat$dead_or_alive,levels=c("Alive","Dead"))
dat <- dat %>% filter(sample_type=='BAL')
    dat <- dat %>% filter(sample_type=='Blood')

#Manually set this based on <21 days | 21 to 42 days | > 42 days 
dat <- dat %>% mutate(is_significant_part1=case_when(
    sample_type=='BAL' & AB_name=='Spike' & AB == 'IgG' ~ 'Y',
    sample_type=='BAL' & AB_name=='SpikeNTD' & AB == 'IgA' ~ 'Y',
    sample_type=='BAL' & AB_name=='SpikeNTD' & AB == 'IgG' ~ 'Y',
    sample_type=='BAL' & AB_name=='SpikeRBD' & AB == 'IgA' ~ 'Y',
    sample_type=='BAL' & AB_name=='SpikeRBD' & AB == 'IgM' ~ 'Y',
    sample_type=='BAL' & AB_name=='SpikeRBD' & AB == 'IgG' ~ 'Y',
    sample_type=='BAL' & AB_name=='TT' & AB == 'IgA' ~ 'Y',
    sample_type=='BAL' & AB_name=='LukS' & AB == 'IgM' ~ 'Y',
    TRUE ~ 'N'
))

for(abname in c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid")){
    ggsave(filename=paste0("lineplot_outcome_blood_",abname,".pdf"),
    # ggsave(filename=paste0("lineplot_outcome_bal_",abname,".pdf"),
    ggplot(dat%>%filter(AB_name==abname,sxs_to_sample_days<=55),aes(x=sxs_to_sample_days,y=AB_value_plus1,group=dead_or_alive,color=dead_or_alive)) +
    geom_vline(xintercept=21, linetype='dashed',color="gray", alpha=0.8)+geom_vline(xintercept=42, linetype='dashed',color="gray", alpha=0.8)+
    geom_smooth(method="loess",se=TRUE,alpha=0.2,level=0.8,size=1.2)+
    scale_x_continuous(limits=c(1,55), breaks=seq(10,55,10),expand = expansion(add = .1))+scale_y_log10(labels=trans_format('log10', math_format(10^.x)))+    
    facet_grid(.~AB)+
    ggtitle(abname)+labs(x="Symptom onset to sample collection (days)", y='Log immunoglobulin level')+
    geom_point(x=12,y=4.7,aes(color=is_significant_part1,shape=is_significant_part1),size=10)+
    scale_color_manual(labels=c("Alive","Dead","Y","N"), limits=c("Alive","Dead","Y","N"), values=c(Alive="#776bcd",Dead="firebrick1",Y="black",N="white"))+scale_fill_manual(values=pal_deadoralive,name="")+scale_shape_manual(labels=c("Y","N"), limits=c("Y","N"), values=c(Y=42,N=32))+
    theme_pubr()+
    theme(legend.position = "none")+theme(axis.line=element_line(size=0.2),plot.title=element_text(hjust=0.5,vjust=1.5,size=txt_size,face="bold"))+theme(axis.title.y=element_text(size=txt_size,vjust=1.5),axis.text.y=element_text(size=txt_size,color="black"))+
    theme(axis.title.x=element_text(vjust=-1.3,size=txt_size),axis.text.x = element_text(size=txt_size,vjust=1,color="black"))+
    theme(strip.background.x=element_blank(),strip.text.x=element_text(color="black",vjust=0.8,size=txt_size))+theme(plot.margin=margin(8,6,12,6,"pt")),
    height=5,width=18)
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Create figure: VIRAL LOAD ANALYSIS
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


updatePearsonCorrelationDf <- function(cur_df,new_row){
    if(nrow(cur_df) == 0){
        cur_df[1,] <- new_row
    }else{
        cur_df <- rbind(cur_df,new_row)    
    } 
    return(cur_df)
}
getPearsonRow <- function(file_touse,timepoint,it,e) {
    tmp_df <- file_touse %>% filter(AB_name==e,AB==it)
    abs_res <- cor.test(tmp_df$AB_value_plus1, tmp_df$updated_absolute_viralload,method="pearson")
    toprint_abs <- c(estimate = abs_res$estimate,p.value = abs_res$p.value,statistic=abs_res$statistic,confint=abs_res$conf.int[1:2])
    norm_res <- cor.test(tmp_df$AB_value_plus1, tmp_df$updated_normalized_viralload,method="pearson")
    toprint_norm <- c(estimate = norm_res$estimate,p.value = norm_res$p.value,statistic=norm_res$statistic,confint=norm_res$conf.int[1:2])    
    return(c(it,e,timepoint,round(toprint_abs[["p.value"]],4),round(toprint_abs[["estimate.cor"]],4),round(toprint_abs[["statistic.t"]],4),round(toprint_norm[["p.value"]],4),round(toprint_norm[["estimate.cor"]],4),round(toprint_norm[["statistic.t"]],4)))
}
make_correlation_figure <- function(pearson_corr_df,figure_name){
    pearson_corr_df$isotype <- factor(pearson_corr_df$isotype , levels=c('IgA','IgM','IgG'))
    pearson_corr_df$epitope <- factor(pearson_corr_df$epitope , levels=c('Nucleocapsid','SpikeRBD','SpikeNTD','Spike'))

    a<-ggplot()+
        geom_hline(yintercept=0, linetype='dashed',color="gray", alpha=0.5)+       
        geom_segment(data=pearson_corr_df, aes(x=epitope, xend=epitope, y=0, yend=as.numeric(absvl_estimate.cor)), color=ifelse(as.numeric(pearson_corr_df$absvl_p.value) <= 0.05,'#ef9b20','gray70'),size=2) + #22a7f0 (blue), #ea5545 (red),#ffb400 (gold), 776bcd (purple), ef9b20 (orange)
        geom_point(data=pearson_corr_df %>% filter(as.numeric(absvl_p.value) <= 0.05), aes(x=epitope,y=as.numeric(absvl_estimate.cor)),color='#ef9b20', size=3,shape=16)+
        coord_flip()+
        theme_pubr()+
        scale_y_continuous(breaks=seq(-1,1,1), limits=c(-1,1))+   
        facet_grid(.~isotype,scales="free",space="free")+
        ylab("Correlation coefficient")+
        theme( legend.position = "none", panel.grid.minor=element_line(color="gray98"), axis.title.y=element_blank(),strip.background = element_blank(),  strip.text.y = element_blank())+
        theme(axis.text.x=element_text(size=txt_size),axis.text.y=element_text(size=txt_size,margin=margin(r=10)), axis.title.x=element_text(size=txt_size))+
        theme(axis.line=element_line(size=1),axis.text.x=element_text(angle = 45, hjust=1,size=txt_size,color="black"),axis.text.y=element_text(size=txt_size,,color="black"),axis.title.y=element_text(size=txt_size,margin=margin(r=10),vjust=1.5),axis.ticks.x=element_blank())+
        theme(strip.background.x=element_blank(), strip.placement = "outside",panel.spacing=unit(0.5,"lines"))+
        theme(strip.text.x=element_text(size=txt_size,vjust=1.9))#,
         ggsave(filename=paste0(figure_name,".pdf"),        
            ggarrange(a,ncol=1,nrow=1),
        height=4,width=9.5)
}



# All BAL samples should have normalized_vl and absolute_vl that are not NA
dat <- useForLinesOrBoxplots(FALSE,FALSE) #FALSE to include all samples (not just baseline), FALSE to not include controls
dat <- dat %>% filter(AB_name %in% c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid"))

# Goal 1: get all blood samples where we have a BAL VL available
dat_bloodfigure <- dat %>% filter(sample_type=='Blood')

#Need external viral load file because for this dataset there was one blood sample where we have BAL viral load data for that day but we did not include in the immunoglobulin levels
bal_vls <- read.csv("allBAL_withviralload.csv", sep=",")
bal_vls <- bal_vls %>% mutate(updated_absolute_viralload=case_when(
    as.numeric(absolute_vl) == 0 ~ as.numeric(absolute_vl)+0.1,
    TRUE ~ as.numeric(absolute_vl)
))
bal_vls <- bal_vls %>% mutate(updated_normalized_viralload=case_when(
    as.numeric(normalized_vl) == 0 ~ as.numeric(normalized_vl)+.00001,
    TRUE ~ as.numeric(normalized_vl)
))

#This file has all blood samples with added columns for the BAL viral loads from that day
dat_bloodfigure_merged <- merge(dat_bloodfigure%>%select(-c("normalized_vl","absolute_vl")), bal_vls%>%select(-c("yl_sample_id","is_sample_id","internal_sample_id")), by=c("study_id","adm_to_sample_days"))
#Only use samples from first 21 days for figure
dat_bloodfigure_merged <- dat_bloodfigure_merged %>% filter((as.numeric(sxs_to_sample_days) >= 0) & (as.numeric(sxs_to_sample_days) < 21))

pearson_corr_df <- data.frame(matrix(ncol=9,nrow=0))
colnames(pearson_corr_df)=c("isotype","epitope","timepoint","absvl_p.value","absvl_estimate.cor","absvl_statistic.t","normvl_p.value","normvl_estimate.cor","norml_statistic.t")
for(it in c("IgA","IgM","IgG")){
    for(e in c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid")){
        new_row <- getPearsonRow(dat_bloodfigure_merged %>% filter(AB_name==e,AB==it), "0-21 days post sxs",it,e)
        pearson_corr_df <- updatePearsonCorrelationDf(pearson_corr_df,new_row)
    }
}

make_correlation_figure(pearson_corr_df,"correlfigure_BALvl_BloodIg_abs")

# Goal 2: get all BAL samples where we have a BAL VL available
dat_balfigure <- dat %>% filter(sample_type=='BAL')

dat_balfigure <- dat_balfigure %>% mutate(updated_absolute_viralload=case_when(
    as.numeric(absolute_vl) == 0 ~ as.numeric(absolute_vl)+0.1,
    TRUE ~ as.numeric(absolute_vl)
))
dat_balfigure <- dat_balfigure %>% mutate(updated_normalized_viralload=case_when(
    as.numeric(normalized_vl) == 0 ~ as.numeric(normalized_vl)+.00001,
    TRUE ~ as.numeric(normalized_vl)
))

#Only use samples from first 21 days for figure
dat_balfigure <- dat_balfigure %>% filter((as.numeric(sxs_to_sample_days) >= 0) & (as.numeric(sxs_to_sample_days) < 21))

pearson_corr_df <- data.frame(matrix(ncol=9,nrow=0))
colnames(pearson_corr_df)=c("isotype","epitope","timepoint","absvl_p.value","absvl_estimate.cor","absvl_statistic.t","normvl_p.value","normvl_estimate.cor","norml_statistic.t")
for(it in c("IgA","IgM","IgG")){
    for(e in c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid")){
        new_row <- getPearsonRow(dat_balfigure %>% filter(AB_name==e,AB==it), "0-21 days post sxs",it,e)
        pearson_corr_df <- updatePearsonCorrelationDf(pearson_corr_df,new_row)
    }
}

make_correlation_figure(pearson_corr_df,"correlfigure_BALvl_BALIg_abs")


### ### ### ### ### ### ###
### ### ### ### ### ### ###
# Ratio VL to IG levels #
### ### ### ### ### ### ###
### ### ### ### ### ### ###

# Ratio at baseline:
dat <- useForLinesOrBoxplots(TRUE,FALSE) #TRUE to include only baseline samples, FALSE to not include controls
dat <- dat %>% filter(AB_name %in% c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid"))

# BAL VL to BAL IgG Spike, SpikeNTD, SpikeRBD, Nucleocapsid
    dat_balfigure <- dat %>% filter(sample_type=='BAL')
    dat_balfigure <- dat_balfigure %>% mutate(updated_absolute_viralload=case_when(
        as.numeric(absolute_vl) == 0 ~ as.numeric(absolute_vl)+0.1,
        TRUE ~ as.numeric(absolute_vl)
    ))
    dat_balfigure <- dat_balfigure %>% mutate(updated_normalized_viralload=case_when(
        as.numeric(normalized_vl) == 0 ~ as.numeric(normalized_vl)+.00001,
        TRUE ~ as.numeric(normalized_vl)
    ))
    fig_dat<-dat_balfigure
# BAL V to Blood IgG 
    dat_bloodfigure <- dat %>% filter(sample_type=='Blood')
    #Need external viral load file because for this dataset there was one blood sample where we have BAL viral load data for that day but we did not include in the immunoglobulin levels
    bal_vls <- read.csv("/Users/barnec03/Research/COVID/SubgroupAnalysis/Immunoglobulins/data_files/allBAL_withviralload.csv", sep=",")
    bal_vls <- bal_vls %>% mutate(updated_absolute_viralload=case_when(
        as.numeric(absolute_vl) == 0 ~ as.numeric(absolute_vl)+0.1,
        TRUE ~ as.numeric(absolute_vl)
    ))
    bal_vls <- bal_vls %>% mutate(updated_normalized_viralload=case_when(
        as.numeric(normalized_vl) == 0 ~ as.numeric(normalized_vl)+.00001,
        TRUE ~ as.numeric(normalized_vl)
    ))
    #This file has all blood samples with added columns for the BAL viral loads from that day
    dat_bloodfigure_merged <- merge(dat_bloodfigure%>%select(-c("normalized_vl","absolute_vl")), bal_vls%>%select(-c("yl_sample_id","is_sample_id","internal_sample_id")), by=c("study_id","adm_to_sample_days"))
    fig_dat <- dat_bloodfigure_merged

fig_dat$AB_name = factor(fig_dat$AB_name,levels=c("Spike","SpikeNTD","SpikeRBD","Nucleocapsid"))
fig_dat$dead_or_alive = factor(fig_dat$dead_or_alive,levels=c("Alive","Dead"))

fig_dat <- fig_dat %>% mutate(ab_to_vl_ratio=log10(AB_value_plus1)/log10(updated_absolute_viralload))
    a<- ggplot(fig_dat, aes(x=AB_name, y=ab_to_vl_ratio,color=factor(dead_or_alive)))+
        geom_boxplot(aes(fill=dead_or_alive,color=dead_or_alive), width=0.5,position=position_dodge(0.7),outlier.shape = NA,alpha=0.3)+
        scale_color_manual(values=pal_deadoralive,name="")+scale_fill_manual(values=pal_deadoralive,name="")+
        facet_grid(cols=vars(AB),scales="free_x",space="free")+
        labs(y='Log10 antibody level/Log10 absolute viral load')+
        theme_pubr()+theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,1,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_text(angle = 45, hjust=1,size=txt_size,color="black"),axis.text.y=element_text(size=txt_size,,color="black"),axis.title.y=element_text(size=txt_size,margin=margin(r=10),vjust=1.5),axis.ticks.x=element_blank())+
        theme(strip.background.x=element_blank(), strip.placement = "outside",panel.spacing=unit(0.5,"lines"))+
        theme(strip.text.x=element_text(size=txt_size))+
        stat_compare_means(aes(label =..p.signif..),label.x=1.5,vjust=0.6,size=10)

        ggsave(filename=paste0("Blood_log10IgToLog10absVL_allABIG_DA_wstats.pdf"),        
            ggarrange(a,ncol=1,nrow=1),
        height=6,width=10)

    fig_dat <- fig_dat %>% mutate(is_significant_part1=case_when(
        sample_type=='BAL' & AB_name=='SpikeRBD' & AB == 'IgA' ~ 'Y',
        sample_type=='BAL' & AB_name=='Nucleocapsid' & AB == 'IgA' ~ 'Y',
        sample_type=='BAL' & AB_name=='Spike' & AB == 'IgG' ~ 'Y',
        sample_type=='BAL' & AB_name=='SpikeRBD' & AB == 'IgG' ~ 'Y',
        sample_type=='BAL' & AB_name=='SpikeNTD' & AB == 'IgG' ~ 'Y',
        sample_type=='BAL' & AB_name=='Nucleocapsid' & AB == 'IgG' ~ 'Y',
        TRUE ~ 'N'
    ))

    a<- ggplot(fig_dat, aes(x=AB_name, y=ab_to_vl_ratio,color=factor(dead_or_alive)))+
        geom_boxplot(aes(fill=dead_or_alive,color=dead_or_alive), width=0.5,position=position_dodge(0.7),alpha=0.3,outlier.shape = NA)+
        scale_color_manual(labels=c("Alive","Dead","Y","N"), limits=c("Alive","Dead","Y","N"), values=c(Alive="#776bcd",Dead="firebrick1",Y="black",N="white"))+scale_fill_manual(values=pal_deadoralive,name="")+scale_shape_manual(labels=c("Y","N"), limits=c("Y","N"), values=c(Y=42,N=32))+
        facet_grid(cols=vars(AB),scales="free_x",space="free")+
        scale_y_log10()+
        labs(y='Ratio AB level/viral load')+
        theme_pubr()+
        theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),axis.title.x=element_blank())+
        theme(plot.margin=margin(l=5))+
        theme(axis.line=element_line(size=1),axis.text.x=element_text(angle = 45, hjust=1,size=txt_size,color="black"),axis.text.y=element_text(size=txt_size,,color="black"),axis.title.y=element_text(size=txt_size,margin=margin(r=10),vjust=1.5),axis.ticks.x=element_blank())+
        theme(strip.background.x=element_blank(), strip.placement = "outside",panel.spacing=unit(0.5,"lines"))+
        theme(strip.text.x=element_text(size=txt_size,vjust=1.9))

        ggsave(filename=paste0("Blood_log10ratio_log10IgToLog10absVL_allABIG_DA_manualstats.pdf"),        
            ggarrange(a,ncol=1,nrow=1),
        height=4.5,width=7.9)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# CASE CONTROL ANALYSIS
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


dat <- useForLinesOrBoxplots(FALSE,TRUE)
txt_size=22
st="Blood"
a <- ggplot(dat %>% filter(sample_type==st) %>% filter(AB_name=='Spike'), aes(x=protocol, y=AB_value_plus1,color=factor(protocol)))+
    geom_boxplot(aes(fill = protocol,color=protocol), alpha = 0.3,width=0.7)  + geom_point(aes(col = protocol))+
    scale_color_manual(values=pal_protocol,name="")+scale_fill_manual(values=pal_protocol,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("Spike")+labs(y='Log Immunoglobulin level')+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size,margin=margin(r=10)),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=1.5),panel.spacing=unit(0.2,"lines"))+#t,r,b,l
    stat_compare_means(aes(label = ifelse(p>5.e-2,..p.signif..,"")),label.x=1.5,vjust=0.5,size=8)+#will be ns one
    stat_compare_means(aes(label = ifelse(p>5.e-2,"",..p.signif..)),label.x=1.5,vjust=1,size=10) #will be other one
    
b <- ggplot(dat %>% filter(sample_type==st) %>% filter(AB_name=='SpikeNTD'), aes(x=protocol, y=AB_value_plus1,color=factor(protocol)))+
    geom_boxplot(aes(fill = protocol,color=protocol), alpha = 0.3,width=0.7)  + geom_point(aes(col = protocol))+
    scale_color_manual(values=pal_protocol,name="")+scale_fill_manual(values=pal_protocol,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("SpikeNTD")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=1.5),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label = ifelse(p>5.e-2,..p.signif..,"")),label.x=1.5,vjust=0.5,size=8)+#will be ns one
    stat_compare_means(aes(label = ifelse(p>5.e-2,"",..p.signif..)),label.x=1.5,vjust=1,size=10) #will be other one
    
c <- ggplot(dat %>% filter(sample_type==st) %>% filter(AB_name=='SpikeRBD'), aes(x=protocol, y=AB_value_plus1,color=factor(protocol)))+
    geom_boxplot(aes(fill = protocol,color=protocol), alpha = 0.3,width=0.7)  + geom_point(aes(col = protocol))+
    scale_color_manual(values=pal_protocol,name="")+scale_fill_manual(values=pal_protocol,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("SpikeRBD")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=1.5),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label = ifelse(p>5.e-2,..p.signif..,"")),label.x=1.5,vjust=0.5,size=8)+#will be ns one
    stat_compare_means(aes(label = ifelse(p>5.e-2,"",..p.signif..)),label.x=1.5,vjust=1,size=10) #will be other one
    
d <- ggplot(dat %>% filter(sample_type==st) %>% filter(AB_name=='Nucleocapsid'), aes(x=protocol, y=AB_value_plus1,color=factor(protocol)))+
    geom_boxplot(aes(fill = protocol,color=protocol), alpha = 0.3,width=0.7)  + geom_point(aes(col = protocol))+
    scale_color_manual(values=pal_protocol,name="")+scale_fill_manual(values=pal_protocol,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("Nucleocapsid")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=1.5),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label = ifelse(p>5.e-2,..p.signif..,"")),label.x=1.5,vjust=0.5,size=8)+#will be ns one
    stat_compare_means(aes(label = ifelse(p>5.e-2,"",..p.signif..)),label.x=1.5,vjust=1,size=10) #will be other one
    
#Use for BAL
    a <- a + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    b <- b + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    c <- c + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    d <- d + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
#Use for Blood
    a <- a + scale_y_log10(limits=c(1,1000000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    b <- b + scale_y_log10(limits=c(1,1000000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    c <- c + scale_y_log10(limits=c(1,1000000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    d <- d + scale_y_log10(limits=c(1,1000000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4

# pdf("boxplot_bal_casecontrol_covidrelated.pdf",width=14.45, height=4.4)
pdf("/boxplot_blood_casecontrol_covidrelated.pdf",width=14.45, height=4.4)
gridExtra::grid.arrange(
    a,b,c,d,
    ncol=4,
    widths=c(3.95,3.5,3.5,3.5),
    heights=c(4.4))
dev.off()

st="BAL"
a <- ggplot(dat %>% filter(sample_type==st) %>% filter(AB_name=='TT'), aes(x=protocol, y=AB_value_plus1,color=factor(protocol)))+
    geom_boxplot(aes(fill = protocol,color=protocol), alpha = 0.3,width=0.7)  + geom_point(aes(col = protocol))+
    scale_color_manual(values=pal_protocol,name="")+scale_fill_manual(values=pal_protocol,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("TT")+labs(y='Log Immunoglobulin level')+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size,margin=margin(r=10)),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=1.5),panel.spacing=unit(0.2,"lines"))+#t,r,b,l
    stat_compare_means(aes(label = ifelse(p>5.e-2,..p.signif..,"")),label.x=1.5,vjust=0.5,size=8)+#will be ns one
    stat_compare_means(aes(label = ifelse(p>5.e-2,"",..p.signif..)),label.x=1.5,vjust=0.8,size=10) #will be other one
    # Works for Blood    
    # stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10)+
    # stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.2,size=7)

b <- ggplot(dat %>% filter(sample_type==st) %>% filter(AB_name=='LukS'), aes(x=protocol, y=AB_value_plus1,color=factor(protocol)))+
    geom_boxplot(aes(fill = protocol,color=protocol), alpha = 0.3,width=0.7)  + geom_point(aes(col = protocol))+
    scale_color_manual(values=pal_protocol,name="")+scale_fill_manual(values=pal_protocol,name="")+
    facet_wrap(c("AB"), nrow=1,strip.position = "bottom",scales="free_x")+
    ggtitle("LukS")+
    theme_pubr()+
    theme(legend.position="none",plot.title=element_text(size=txt_size,hjust=0.5,vjust=2.5),axis.ticks=element_line(size=0.3),plot.margin=margin(10,6,5.5,5.5,"pt"),axis.title.x=element_blank(),axis.line=element_line(size=0.2),axis.text.x=element_blank(),axis.text.y=element_text(size=txt_size),axis.title.y=element_text(size=txt_size),axis.ticks.x=element_blank())+
    theme(strip.background.x=element_blank(), strip.placement = "outside",strip.text.x=element_text(size=txt_size,vjust=1.5),panel.spacing=unit(0.2,"lines"))+
    theme(axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+        
    stat_compare_means(aes(label = ifelse(p>5.e-2,..p.signif..,"")),label.x=1.5,vjust=0.5,size=8)+#will be ns one
    stat_compare_means(aes(label = ifelse(p>5.e-2,"",..p.signif..)),label.x=1.5,vjust=0.8,size=10) #will be other one
    # Works for Blood    
    # stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,"",..p.signif..)),label.x=1.5,vjust=0.6,size=10)+
    # stat_compare_means(aes(label =ifelse(after_stat(p.format)>0.05,..p.signif..,"")),label.x=1.5,vjust=0.2,size=7)
#For BAL
    a <- a + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    b <- b + scale_y_log10(limits=c(.1,100000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
#For Blood
    a <- a + scale_y_log10(limits=c(1,1000000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4
    b <- b + scale_y_log10(limits=c(1,1000000), labels=trans_format('log10', math_format(10^.x)))#limits=c(5,100000), ) 0 to 4

pdf("boxplot_bal_casecontrol_notcovidrelated.pdf",width=7.45, height=4.4)
# pdf("boxplot_blood_casecontrol_notcovidrelated.pdf",width=7.45, height=4.4)
gridExtra::grid.arrange(
    a,b,
    ncol=2,
    widths=c(3.95,3.5),
    heights=c(4.4))
dev.off()





