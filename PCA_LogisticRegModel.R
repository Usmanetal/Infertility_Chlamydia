




library(readxl)
library(janitor)
library(tidyverse)
library(finalfit)
library(dplyr)
library(MASS)
library(factoextra)
library(rstatix)
library(ggpubr)
Chlamydia_Venn <- read_excel("C:/Users/Usman Ola/Downloads/Chlamydia Trichomatis.xlsx")
Chlamydia_logistic_reg <- clean_names(Chlamydia_Venn)
names(Chlamydia_logistic_reg) # All variables names are automagically cleaned!
newdata<-Chlamydia_logistic_reg[!Chlamydia_logistic_reg$how_long_7=="Nil",]%>%
  mutate(age_cat=cut(age,
                     breaks = c(0,35,Inf),
                     labels = c("<35",">35")),
         Marital_dur=ifelse(how_long_7=="1","First_Child","Children"),
         Marital_cat=cut(as.numeric(how_long_7),
                         c(0,20,Inf),
                         c("≤20",">20")),
         Times_mar=case_when(how_times=="Fouth"~"second",
                             how_times=="Third"~"second",
                             how_times=="First"~"first",
                             how_times=="Second"~"second"),
         wom_edu=str_replace(level_edu,"p","P"),
         husband_occ=case_when(
           husb_occ=="Leccturer"~"Lecturer",
           husb_occ=="Teacher"~"Lecturer",
           husb_occ=="Civil Serv"~"Civil_ser",
           husb_occ=="Lecturer"~"Lecturer",
           husb_occ=="Civil serv"~"Civil_ser",
           husb_occ=="Publicserv"~"Civil_ser",
           husb_occ=="Civilserv"~"Civil_ser",
           husb_occ=="Bussiness"~"Bussines",
           husb_occ=="Business"~"Bussines",
           husb_occ=="Carpenter"~"Artisan",
           husb_occ=="Artisan"~"Artisan",
           husb_occ=="Farmer"~"Farmer",
           TRUE~"others"),
         fam_type=str_replace_all(family_typ,c("mono"="Mono")),
         woman_occ=case_when(
           occup=="Teacher"~"Teacher",
           occup=="teacher"~"Teacher",
           occup=="Civil Serv"~"Civil_ser",
           occup=="House wife"~"House_wife",
           occup=="Civil serv"~"Civil_ser", 
           occup=="Civilserv"~"Civil_ser",
           occup=="Student"~"Student",
           occup=="Tailor"~"Artisan",
           occup=="Artisan"~"Artisan", 
           occup=="Farmer"~"Farmer",
           TRUE~"others"),woman_occ1=ifelse(woman_occ=="Teacher","Teacher","other"),
         woman_occ2=str_replace_all(woman_occ,c("Farmer"="others","Student"="others","Teacher"="others","Artisan"="others")),
         woman_occ2f=factor(woman_occ2,levels=c("others","House_wife","Civil_ser")),
         any_preg=recode(any_preg,"No"="Primary_infert",
                         "Yes"="secondary_infert"),
         how_many=recode(
           how_many,
           "Eight"=8,
           "Five"=5,
           "Four"=4,
           "Nil"=0,
           "Once"=1,
           "One"=1,
           "Seven"=7,
           "Six"=6,
           "Ten"=10,
           "Three"=3,
           "Thrice"=3,
           "Twice"=2,
           "Twince"=2,
           "Two"=2,
           "10"=10,
           "2"=2,
           "3"=3),
         Pregnancy=cut(how_many,c(-Inf,0.9,1,5,10),
                       c("No_preg","once","2-5","6-10")), 
         to_term=recode(
           to_term,
           "Eight"=8,
           "Five"=5,
           "Four"=4,
           "Nil"=0,
           "Once"=1,
           "One"=1,"0ne"=1,"one"=1,
           "Seven"=7,
           "Six"=6,
           "Ten"=10,
           "Three"=3,
           "Thrice"=3,
           "Twice"=2,
           "Twince"=2,
           "Two"=2,
           "10"=10,
           "2"=2,
           "3"=3,"Nine"=9),
         term_cat=cut(to_term,c(-Inf,0.9,1,5,10),right=TRUE,
                      c("no_child","One_child","2-5","6-10"),.drop=FALSE) ,
         age_of_last=as.numeric(str_replace(age_of_last,"Nil","0")),
         ageoflast_cat=cut(age_of_last,c(-Inf,0.9,1,5,10,15,20),
                           c("none","One_year","1-5","6-10","11-15","16-20")),
         family_pla,
         which_one,
         how_long_20=recode(
           how_long_20,"1"=1,"2"=2,"3"=3,"5"=5,"6"=6,"Nil"=0 ,"<6M"=1,"NiL"=0),
         Vaginal_discharge=str_replace_all(color,c("NIL"="Nil","Brownish"="Reddish")),
         vaginal_discharge1=str_replace_all(Vaginal_discharge,c("Brownish"="Yellowish","Reddish"="Yellowish")),
         il_10_quatile=ntile(il_10_avrg,3),
         il_10_quatile_cat=recode(il_10_quatile,"1"="firstq","2"="sec","3"="thir"),
         ifn_gam_quatile=ntile(ifn_gamm_avrg,3),
         ifn_gam_quatile_cat=recode(ifn_gam_quatile,"1"="firstq","2"="sec","3"="thir"),
         hsp_60_quatile=ntile(hsp_60_avrg,3),chlamydia_pos=ifelse(ig_g=="Pos"|ig_m=="Pos","Pos","Neg"),
         hsp_60_quatile_cat=recode(hsp_60_quatile,"1"="firstq","2"="sec","3"="thir")
         
  )


## **R Markdown of cytokine profile in infertile women: focus on dimension \n reduction**

### **Principal Component Analysis (PCA) for dimension reduction**
### **PCA for cluster recognition**
## **Preparing data for dimension and variable contribution**
### **Scree Plot**



cytokine_response<-newdata[newdata$group=="Infertile",c("il_10_avrg","ifn_gamm_avrg","hsp_60_avrg")]%>%rename("IL-10"="il_10_avrg","paste(IFN-gamma)"="ifn_gamm_avrg","Ct-Hsp60"="hsp_60_avrg")

# Run PCA
cytokine_pca<- prcomp(cytokine_response,scale=TRUE)

# Scree plot of variance

fviz_eig(cytokine_pca,
         addlabels = TRUE,
         ylim=c(0,70),
         main = "Cytokine screePlot"
)



### PCA Plot


# Biplot with Default Settings
fviz_pca_biplot(cytokine_pca,parse=TRUE)

# Biplot with labeled variables
fviz_pca_biplot(cytokine_pca,
                label = "var", parse=TRUE)


### Pattern Recognition


# Biplot with customized colored groups and variables
Ct_pos_neg_pca<-fviz_pca_biplot(cytokine_pca,parse=TRUE,repel = TRUE,
                                label="var",addEllipses = TRUE,title="Cytokine PCA: Infertility",
                                habillage=newdata[newdata$group=="Infertile",]$chlamydia_pos  ,col.var="black")+labs(x="PCA2 (57.1%)",y="PCA2 (29.8%)")+theme(legend.title  = element_blank())+scale_fill_manual(labels=c("Ct-negative","Ct-positive"),
                                                                                                                                                                                                                 values=c("orange","purple"))+scale_shape_manual(labels=c("Ct-negative",
                                                                                                                                                                                                                                                                          "Ct-positive"),values=c(19,17))+scale_color_manual(labels=c("Ct-negative",    "Ct-positive"),values=c("orange","purple"))


## Principal Component Analysis for cytokines in Fertile and Infertile Women

# Cytokine Response in fertile and infertile women
cytokine_response_g<-newdata[,c("il_10_avrg","ifn_gamm_avrg","hsp_60_avrg")] %>%rename("IL-10"="il_10_avrg","paste(IFN-gamma)"="ifn_gamm_avrg","Ct-Hsp60"="hsp_60_avrg")

# Run PCA
cytokine_pca_g<- prcomp(cytokine_response_g,scale=TRUE)

# Scree plot of variance

fviz_eig(cytokine_pca_g,
         addlabels = TRUE,
         ylim=c(0,70),
         main = "Cytokine screePlot"
)

# Biplot with Default Settings
fviz_pca_biplot(cytokine_pca_g)

# Biplot with labeled variables
fviz_pca_biplot(cytokine_pca_g,
                label = "var")

# Biplot with customized colored groups and variables
fviz_pca_biplot(cytokine_pca_g,
                label = "var",
                habillage = newdata$chlamydia_pos)

# Biplot with customized colored groups and variables
Ct_Pri_sec_both<-fviz_pca_biplot(cytokine_pca_g,
                                 label="var",addEllipses = TRUE,title="Cytokine PCA: both Fertile and Infertile",
                                 habillage=newdata$chlamydia_pos,parse=TRUE,
                                 col.var="black")+scale_color_manual(values=c("orange","purple"))+labs(x="PCA1 (56.3%)",y="PCA2 (29.9%)")+theme(legend.title  = element_blank())+scale_fill_manual(labels=c("Ct-negative","Ct-positive"),
                                                                                                                                                                                                   values=c("orange","purple"))+scale_shape_manual(labels=c("Ct-negative",
                                                                                                                                                                                                                                                            "Ct-positive"),values=c(19,17))+scale_color_manual(labels=c("Ct-negative",    "Ct-positive"),values=c("orange","purple"))

Figure_pca<-ggarrange(Ct_pos_neg_pca,Ct_Pri_sec_both,labels = c("A",
                                                                "B"),common.legend = TRUE,legend = "bottom")


## **Cytokine levels in Fertile and Infertile subset by their Chlamydial positivity**


pwc <- newdata %>%
  #group_by(chlamydia_pos)%>%
  dunn_test(il_10_avrg~chlamydia_pos,p.adjust.method = "none")
pwc <- pwc %>% add_y_position()
pwc
ggplot(data=newdata,aes(x=any_preg,y=hsp_60_avrg))+
  geom_violin(trim = FALSE,col="blue",size=1.5)+
  facet_wrap(~group,nrow=1)

Interleukin_10<-ggplot(data=newdata,aes(x=group,y=il_10_avrg))+
  geom_violin(trim = FALSE,col="blue",size=1.5)+
  geom_boxplot(width=0.1)+
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  stat_pvalue_manual(pwc, tip.length = 0.02,hide.ns = FALSE,
                     bracket.size = 1,size=7,label = "p.adj.signif",color="darkblue")+facet_wrap(~chlamydia_pos)

pwc1 <- newdata %>%group_by(chlamydia_pos)%>%
  dunn_test(ifn_gamm_avrg~group,p.adjust.method = "none")
pwc1 <- pwc1 %>% add_y_position()
pwc1
Interferon_alpha<-ggplot(data=newdata,aes(x=group,y=ifn_gamm_avrg))+
  geom_violin(trim = FALSE,col="blue",size=1.5)+
  geom_boxplot(width=0.1)+facet_wrap(~chlamydia_pos)+
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  stat_pvalue_manual(pwc1, tip.length = 0.02,hide.ns = FALSE,
                     bracket.size = 1,size=7,label = "p.adj.signif",color="darkblue" )

pwc2 <- newdata %>%group_by(chlamydia_pos)%>%
  dunn_test(hsp_60_avrg~group,p.adjust.method = "none")
pwc2 <- pwc2 %>% add_y_position()
pwc2
Heat_shock_Prot<-newdata%>%group_by(group)%>%ggplot(aes(x=group,y=hsp_60_avrg))+
  geom_violin(trim = FALSE,col="blue",size=1.5)+
  geom_boxplot(width=0.1)+facet_wrap(~chlamydia_pos)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  stat_pvalue_manual(pwc2, tip.length = 0.02,hide.ns = TRUE,
                     bracket.size = 1,size=7,label = "p.adj.signif",color="darkblue" )

# Cytokine Levels of Fertile and Infertile women 
ggarrange(Interleukin_10,Interferon_alpha,Heat_shock_Prot,nrow = 1)



### Cytokine levels Infertility data subset by Primary or Secondary Infertility and their Chlamydial sensitivity


pwc3 <- newdata[newdata$group=="Infertile",] %>%
  group_by(chlamydia_pos)%>%
  dunn_test(il_10_avrg~any_preg,p.adjust.method = "none")
pwc3 <- pwc3 %>% add_y_position()
pwc3
Interleukin<-newdata[newdata$group=="Infertile",]%>%
  ggplot(aes(x=any_preg,y=il_10_avrg))+
  geom_violin(trim = FALSE,col="blue",size=1.5)+
  geom_boxplot(width=0.1)+facet_wrap(~chlamydia_pos,nrow=1)
# geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
stat_pvalue_manual(pwc3, tip.length = 0.02,hide.ns = FALSE,
                   bracket.size = 1,size=7,label = "p.adj.signif",color="darkblue" )

pwc5 <- newdata[newdata$group=="Infertile",] %>%
  group_by(chlamydia_pos)%>%
  dunn_test(ifn_gamm_avrg~any_preg,p.adjust.method = "none")
pwc5 <- pwc5 %>% add_y_position()
pwc5
Interferon<-newdata[newdata$group=="Infertile",]%>%
  ggplot(aes(x=any_preg,y=ifn_gamm_avrg))+
  geom_violin(trim = FALSE,col="blue",size=1.5)+
  geom_boxplot(width=0.1)+facet_wrap(~chlamydia_pos,nrow=1)
# geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
stat_pvalue_manual(pwc5, tip.length = 0.02,hide.ns = FALSE,
                   bracket.size = 1,size=7,label = "p.adj.signif",color="darkblue" )

pwc6 <- newdata[newdata$group=="Infertile",] %>%
  group_by(chlamydia_pos)%>%
  dunn_test(hsp_60_avrg~any_preg,p.adjust.method = "none")
pwc6 <- pwc6 %>% add_y_position()
pwc6
Hsp<-newdata[newdata$group=="Infertile",]%>%
  ggplot(aes(x=any_preg,y=hsp_60_avrg))+
  geom_violin(trim = FALSE,col="blue",size=1.5)+
  geom_boxplot(width=0.1)+facet_wrap(~chlamydia_pos,nrow=1)
# geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
stat_pvalue_manual(pwc6, tip.length = 0.02,hide.ns = FALSE,
                   bracket.size = 1,size=7,label = "p.adj.signif",color="darkblue" )

ggarrange(Interleukin,Interferon,Hsp,nrow = 1)


## Socio_Demographic Factor of Fertile in Infertile women of our study population


age_cat<-Chlamydia_logistic_reg%>%group_by(group)%>%
  mutate(age1=cut(age,
                  breaks = c(0,35,Inf),
                  labels = c("<35",">35")))%>%
  count(age1)%>% 
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="age1")
#%>%chisq.test(correct=FALSE)

#R code 1
Marital_Duration_couples_mmm<-Chlamydia_logistic_reg[!Chlamydia_logistic_reg$how_long_7=="Nil",]%>%
  reframe(median=median(as.numeric(how_long_7)),
          median=median(as.numeric(how_long_7)),mode=mode(as.numeric(how_long_7)),
          qs=list(quantile(as.numeric(how_long_7))))%>%
  unnest_wider(qs)

Marital_Duration_couples_mmm<-Chlamydia_logistic_reg%>%
  mutate(how_long_7=ifelse(how_long_7=="1","First_Child","Children"))%>%
  group_by(group)%>%
  count(how_long_7)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="how_long_7")

#R code 2
#Chlamydia_logistic_reg[-100,Chlamydia_logistic_reg$how_long_7]%>%
# select((as.numeric(Chlamydia_logistic_reg$how_long_7)))%>%
#summarise_at(
# vars(.),list(min=min,max=max,
#             median=median,Q1=quantile(.,0.25),Q3=quantile(.,0.75)))
# the above code works only when cleaned-up as numeric

Marital_Duration_Couples_Categ<-Chlamydia_logistic_reg[!Chlamydia_logistic_reg$how_long_7=="Nil",]%>%
  group_by(group)%>%
  mutate(categ=cut(as.numeric(how_long_7),c(0,10,20,Inf),c("<10","11-20",">20")))%>%
  count(categ)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="categ")


#no. of life time sex partner
no_times_married<-Chlamydia_logistic_reg[!Chlamydia_logistic_reg$how_times=="Nil",]%>%
  group_by(group)%>%
  mutate(Times_mar=case_when(how_times=="Fouth"~"second",how_times=="Third"~"second",how_times=="First"~"first",how_times=="Second"~"second"))%>%
  count(Times_mar)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="Times_mar")

woman_Education<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  mutate(wom_edu=str_replace(level_edu,"p","P"))%>%count(wom_edu)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="wom_edu")

Husban_Occupation<-Chlamydia_logistic_reg%>%group_by(group)%>%
  mutate(husband_occ=case_when(
    husb_occ=="Leccturer"~"Lecturer",
    husb_occ=="Teacher"~"Lecturer",
    husb_occ=="Civil Serv"~"Civil_ser",
    husb_occ=="Lecturer"~"Lecturer",
    husb_occ=="Civil serv"~"Civil_ser",
    husb_occ=="Publicserv"~"Civil_ser",
    husb_occ=="Civilserv"~"Civil_ser",
    husb_occ=="Bussiness"~"Bussines",
    husb_occ=="Business"~"Bussines",
    husb_occ=="Carpenter"~"Artisan",
    husb_occ=="Artisan"~"Artisan",
    husb_occ=="Farmer"~"Farmer",
    TRUE~"others"))%>%
  count(husband_occ)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="husband_occ")

Fam_type<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  mutate(fam_type=str_replace_all(family_typ,c("mono"="Mono")))%>%
  count(fam_type)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="fam_type")

Woman_occupation<-Chlamydia_logistic_reg%>%group_by(group)%>%
  mutate(woman_occ=case_when(
    occup=="Teacher"~"Teacher",
    occup=="teacher"~"Teacher",
    occup=="Civil Serv"~"Civil_ser",
    occup=="House wife"~"House_wife",
    occup=="Civil serv"~"Civil_ser", 
    occup=="Civilserv"~"Civil_ser",
    occup=="Student"~"Student",
    occup=="Tailor"~"Artisan",
    occup=="Artisan"~"Artisan", 
    occup=="Farmer"~"Farmer",
    TRUE~"others"))%>%
  count(woman_occ)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="woman_occ")%>%mutate(Fertile=replace_na(Fertile,0))

Infertility_cat<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  count(any_preg)%>%
  mutate(any_preg=recode(any_preg,"No"="Primary_infert",
                         "Yes"="secondary_infert"))%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="any_preg")



History_pregnancy<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  mutate(how_many=recode(
    how_many,
    "Eight"=8,
    "Five"=5,
    "Four"=4,
    "Nil"=0, # changed from 0 to 12 to allow for analysis
    "Once"=1,
    "One"=1,
    "Seven"=7,
    "Six"=6,
    "Ten"=10,
    "Three"=3,
    "Thrice"=3,
    "Twice"=2,
    "Twince"=2,
    "Two"=2,
    "10"=10,
    "2"=2,
    "3"=3),
    grp=cut(how_many,c(-Inf,0.9,1,5.0,10),
            c("no_Child","one_child","1-5","6-10")))%>%count(grp)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="grp")%>%mutate(Fertile=replace_na(Fertile,0))
# the 'Nil' has been changed to for 0 to 12 to allow for chi-sq.test analysis

To_term<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  mutate(to_term=recode(
    to_term,
    "Eight"=8,
    "Five"=5,
    "Four"=4,
    "Nil"=0,
    "Once"=1,
    "One"=1,"0ne"=1,"one"=1,
    "Seven"=7,
    "Six"=6,
    "Ten"=10,
    "Three"=3,
    "Thrice"=3,
    "Twice"=2,
    "Twince"=2,
    "Two"=2,
    "10"=10,
    "2"=2,
    "3"=3,"Nine"=9),
    grp=cut(to_term,c(-Inf,0.9,1,5,10),right = TRUE,ordered_result = TRUE,
            c("no_Child","one_child","1-5","6-10")))%>%count(grp)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="grp")%>%mutate(Fertile=replace_na(Fertile,0))
#%>%chisq.test()

Last_age_of_child<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  mutate(age_of_last=as.numeric(str_replace(age_of_last,"Nil","0")),
         grp=cut(age_of_last,c(-Inf,0.9,1,5,10,15,20),
                 c("none","One_year","1-5","6-10","11-15","16-20"),.drop=FALSE))%>%
  count(grp)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="grp")%>%
  mutate(Fertile=replace_na(Fertile,0),Infertile=replace_na(Infertile,0))
#%>%chisq.test()

Family_plan<-Chlamydia_logistic_reg%>%
  group_by(group)%>%count(family_pla)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="family_pla")%>%
  mutate(Fertile=replace_na(Fertile,0),
         Infertile=replace_na(Infertile,0))

What_type_of_contraceptive<-Chlamydia_logistic_reg%>%
  group_by(group)%>%count(which_one)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="which_one")%>%
  mutate(Fertile=replace_na(Fertile,0),
         Infertile=replace_na(Infertile,0))

How_long_use_contra<-Chlamydia_logistic_reg%>%group_by(group)%>%
  mutate(how_long_20=recode(how_long_20,"1"=1,"2"=2,"3"=3,"5"=5,"6"=6,"Nil"=0 ,"<6M"=1,"NiL"=0))%>%
  mutate(how_long_20=cut(how_long_20,c(-Inf,0.9,1,6),c("none","one_year","2-6")))%>%
  count(how_long_20)%>%pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="how_long_20")

#Vaginal disc or menstrual Disc/cup originally used to prevent menstrual bleeding
#as an Eco-friendly alternative to tampon is increasing used for fertility
#(prevention of sperm backflow from the vagina)
# as a second-use
Vaginal_disc<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  count(vaginal_disc)%>%
  pivot_wider(names_from=group,values_from=n)%>% remove_rownames %>%
  column_to_rownames(var="vaginal_disc") 

Vaginal_discharge<-Chlamydia_logistic_reg%>%group_by(group)%>%
  count(color=str_replace_all(color,c("NIL"="Nil")))%>%
  pivot_wider(names_from=group,values_from=n)%>% remove_rownames %>%
  column_to_rownames(var="color") 

Odour<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  count(odour)%>%
  pivot_wider(names_from=group,values_from=n)%>% remove_rownames %>%
  column_to_rownames(var="odour")

#Pelvic Pain and discomfort
Discomfort<-Chlamydia_logistic_reg%>%
  group_by(group)%>%
  count(discomfort)%>%
  pivot_wider(names_from=group,values_from=n)%>% remove_rownames %>%
  column_to_rownames(var="discomfort")%>%mutate(Infertile=replace_na(Infertile,0))

## Quatiles categorization of Interleukin 10

Inter_10<-Chlamydia_logistic_reg%>%group_by(group)%>%
  mutate(il_10_quatile=ntile(il_10_avrg,3))%>%
  count(il_10_quatile)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="il_10_quatile")# %>%
# select(Infertile,Fertile) # alternative to "column_to_rownames"
#chisq.test()
# xtabs(~Infertile+Fertile,data=.) # create a matrix table

## Quatiles to categorize Interferon Alpha

Inter_Alpha<-Chlamydia_logistic_reg%>%group_by(group)%>%
  mutate(ifn_gamm_avrg=ntile(ifn_gamm_avrg,3))%>%
  count(ifn_gamm_avrg)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="ifn_gamm_avrg")#%>%
# select(Infertile,Fertile) # alternative to "column_to_rownames"
#chisq.test()
# xtabs(~Infertile+Fertile,data=.) # create a matrix table

## Categorization of cHSP-60 (Chlamydia Heat Shock Protein)

Chsp_60<-Chlamydia_logistic_reg%>%group_by(group)%>%
  mutate(hsp_60_avrg=ntile(hsp_60_avrg,3))%>%
  count(hsp_60_avrg)%>%
  pivot_wider(names_from=group,values_from=n)%>%
  column_to_rownames(var="hsp_60_avrg")#%>%
# select(Infertile,Fertile) # alternative to "column_to_rownames"
#chisq.test()
# xtabs(~Infertile+Fertile,data=.) # create a matrix table

# Chlamydia anitgenemia

Chlamydia_logistic_reg%>%
  group_by(group)%>%
  mutate(new=ifelse(ig_g=="Pos"|ig_m=="Pos","Pos","Neg"))%>%count(new)

summary(Chlamydia_logistic_reg$il_10_avrg)
summary(Chlamydia_logistic_reg$ifn_gamm_avrg)
summary(Chlamydia_logistic_reg$hsp_60_avrg) 

# Final build up of the Socio-demographic factors surrounding reproductive pathology

Socio_eco<-bind_rows(list(
  age=age_cat,Marital=Marital_Duration_Couples_Categ,
  no_times_mar=no_times_married,
  wom_Edu=woman_Education,
  Husban_Occupation=Husban_Occupation,
  Fam_type=Fam_type,
  Woman_occupation=Woman_occupation,
  Infertility_cat=Infertility_cat,
  History_pregnancy=History_pregnancy,
  To_term=To_term,
  Last_age_of_child=Last_age_of_child,
  Family_plan=Family_plan,
  What_type_of_contraceptive=What_type_of_contraceptive,
  How_long_use_contra=How_long_use_contra,
  Vaginal_disc=Vaginal_disc,
  Vaginal_discharge=Vaginal_discharge,
  Odour=Odour,
  Discomfort=Discomfort,
  Inter_10=Inter_10,
  Inter_Alpha=Inter_Alpha,
  Chsp_60=Chsp_60),.id="Variable")%>%
  rownames_to_column()

# Final sociodemographic data of our study group

Socio_eco%>%mutate(P_Value=recode
                   (Variable,
                     "age"="0.09419",
                     "Marital"="0.0078736",
                     "no_times_mar"="0.0005033",
                     "wom_Edu"="0.009383",
                     "Husban_Occupation"="0.07422",
                     "Fam_type"="0.08756",
                     "Woman_occupation"="<0.0001",
                     "Infertility_cat"="-",
                     "History_pregnancy"="<0.0001",
                     "To_term"="<0.0001",
                     "Last_age_of_child"="=0.001",
                     "Family_plan"="<0.0001",
                     "What_type_of_contraceptive"="0.0378",
                     "How_long_use_contra"="0.02004",
                     "Vaginal_disc"="<0.0001",
                     "Vaginal_discharge"="0.0001",
                     "Odour"="<0.0001","Discomfort"="<0.0001",
                     "Inter_10"="n.s","Inter_Alpha"="n.s","Chsp_60"="n.s"))%>%
  mutate(rowname=str_replace_all(rowname,c("...69"="Quatile1",
                                           "...70"="Quatile2",
                                           "...71"="Quatile3",
                                           "...72"="Quatile_1",
                                           "...73"="Quatile_2",
                                           "...74"="Quatile_3",
                                           "...75"="Quatile.1",
                                           "...76"="Quatile.2",
                                           "...77"="Quatile.3"))
  )


## Exploratory Logistic Data Analysis
### Infertiliy
### Primary and Secondary Infertility


knitr::opts_chunk$set(echo = TRUE)
# Logistic regression Plot for model1

new_logistic<-newdata%>%
  rename("IgM"="ig_m","IgG"="ig_g","Vaginal.discharge.positive"="vaginal_disc",
         "Odour"="odour","Marital.category"="Marital_cat",
         "Number.of.times.maried"="Times_mar","Woman.occupation"="woman_occ1",
         "Age.category"="age_cat","Woman.level.education"="wom_edu",
         "number.of.children"="how_many","Vaginal.discharge.color"="Vaginal_discharge","History.of.Family.Planning"="family_pla",
         "Age"="age",
         "Ct.Heat.shock.protien"="hsp_60_quatile_cat",
         "Family.type"="fam_type",
         "Ct.Positive"="chlamydia_pos")%>%mutate(Ct.Heat.shock.protien=recode(Ct.Heat.shock.protien,"firstq"="Q1","sec"="Q2","thir" ="Q3"))
new_logistic$group<-as.factor(new_logistic$group)


explanatory=c("IgM","IgG","Vaginal.discharge.positive","Odour","Marital.category","Number.of.times.maried","Woman.occupation","Age.category","Woman.level.education","number.of.children")
#%>%mutate(explanatory=str_replace(explanatory,".",""))

dependent = "group"

Model1<-new_logistic %>%mutate(Woman.level.education=recode(Woman.level.education,"N"="None","P"="Primary","S"="Secondary","T"="Tertiary"),Woman.occupation=factor(Woman.occupation,levels = c("other","Teacher")))%>%
  or_plot(dependent, explanatory, table_text_size=3, 
          title_text_size=12,breaks = c(-1,0,1),dependent_label = "Predictors of Infertility (Model1)",remove_ref = TRUE,
          plot_opts=list(xlab("Beta, 95% CI"),
                         theme(axis.title = element_text(size=12))))

newdata$group<-as.factor(newdata$group)
# explanatory = c( "ig_m",
#                    "ig_g",
#                    "vaginal_disc",
#                    "odour",
#                    "Marital_cat",
#                    "Times_mar",
#                    "woman_occ1","age_cat","wom_edu",
#                    "how_many"
# )
dependent = "group"
# newdata %>%
#   or_plot(dependent, explanatory, table_text_size=3, 
#           title_text_size=12,
#           plot_opts=list(xlab("Beta, 95% CI"),
#                          theme(axis.title = element_text(size=12))))

# Logistic regression Plot for model2 (including Vaginal_discharge and removing Vaginal disc)
explanatory=c("IgM","IgG","Vaginal.discharge.color","Odour","Marital.category","Number.of.times.maried","Woman.occupation","Age.category","Woman.level.education","number.of.children")

Model2<-new_logistic %>%mutate(Woman.level.education=recode(Woman.level.education,"N"="None","P"="Primary","S"="Secondary","T"="Tertiary"),Woman.occupation=factor(Woman.occupation,levels = c("other","Teacher")))%>%
  or_plot(dependent, explanatory, table_text_size=3, 
          title_text_size=12,breaks = c(-1,0,1),dependent_label = "Predictors of Infertility (Model2)",remove_ref = TRUE,
          plot_opts=list(xlab("Beta, 95% CI"),
                         theme(axis.title = element_text(size=12))))

figure_Model<-ggarrange(Model1,Model2,nrow = 2,ncol = 1)

# explanatory1 = c( "ig_m",
#                  "ig_g",
#                  "Vaginal_discharge",
#                  "odour",
#                  "Marital_cat",
#                  "Times_mar",
#                  "woman_occ1","age_cat","wom_edu",
#                  "how_many"
# )
# 
# dependent1 = "group"
# newdata %>%
#   or_plot(dependent, explanatory, table_text_size=3, 
#           title_text_size=12,
#           plot_opts=list(xlab("Beta, 95% CI"), 
#                          theme(axis.title = element_text(size=12))))

# Socio-demographic risk factors for Infertility type in the study population


# Model2<-new_logistic %>%mutate(Woman.level.education=recode(Woman.level.education,"N"="None","P"="Primary","S"="Secondary","T"="Tertiary"),Woman.occupation=factor(Woman.occupation,levels = c("other","Teacher")))%>%
#     or_plot(dependent, explanatory, table_text_size=3, 
#             title_text_size=12,breaks = c(-1,0,1),dependent_label = "Predictors of Infertility (Model2)",remove_ref = TRUE,
#             plot_opts=list(xlab("Beta, 95% CI"),
#                            theme(axis.title = element_text(size=12))))

newdata$any_preg<-factor(newdata$any_preg,levels=c("Primary_infert","secondary_infert"))

newdata[newdata$group=="Infertile",]%>%
  glm(factor(any_preg)~family_pla+
        age+hsp_60_quatile_cat+
        fam_type+odour+
        chlamydia_pos,
      data = .,family="binomial")%>%
  summary()
# model to predict Secondary infertility
explanatory = c( "family_pla",
                 "age",
                 "hsp_60_quatile_cat",
                 "fam_type",
                 "odour",
                 "chlamydia_pos"
)

explanatory=c("History.of.Family.Planning",
              "Age","Ct.Heat.shock.protien",
              "Family.type",
              "Odour",
              "Ct.Positive")
dependent = "any_preg"

new_logistic$any_preg<-factor(new_logistic$any_preg,levels=c("Primary_infert","secondary_infert"))

# newdata[newdata$group=="Infertile",] %>%
#   or_plot(dependent, explanatory, table_text_size=3, 
#           title_text_size=12,
#           plot_opts=list(xlab("Beta, 95% CI"), 
#                          theme(axis.title = element_text(size=12))))

Infertility_Model1<-new_logistic[new_logistic$group=="Infertile",] %>%
  or_plot(dependent, explanatory, table_text_size=3, 
          title_text_size=12,breaks = c(-1,0,1),dependent_label = "Predictors of Secondary Infertility",remove_ref = TRUE,
          plot_opts=list(xlab("Beta, 95% CI"),
                         theme(axis.title = element_text(size=12))))

# Model (sociodemographic) to predict Primary Infertility

# model to predict Secondary infertility
explanatory = c( "family_pla",
                 "age",
                 "hsp_60_quatile_cat",
                 "fam_type",
                 "odour",
                 "chlamydia_pos"
)

explanatory=c("History.of.Family.Planning",
              "Age","Ct.Heat.shock.protien",
              "Family.type",
              "Odour",
              "Ct.Positive")
dependent = "any_preg"

newdata$any_preg<-factor(newdata$any_preg,levels=c("secondary_infert","Primary_infert"))
dependent = "any_preg"

new_logistic$any_preg<-factor(new_logistic$any_preg,levels=c("secondary_infert","Primary_infert"))
dependent = "any_preg"

Fertility_Model2<-new_logistic[new_logistic$group=="Infertile",] %>%
  or_plot(dependent, explanatory, table_text_size=3, 
          title_text_size=12,breaks = c(-1,0,1),dependent_label = "Predictors of Primary Infertility",remove_ref = TRUE,
          plot_opts=list(xlab("Beta, 95% CI"),
                         theme(axis.title = element_text(size=12))))

figure_Model2<-ggarrange(Fertility_Model2,Infertility_Model1,ncol = 1,nrow = 2)

# newdata[newdata$group=="Infertile",] %>%
#   or_plot(dependent, explanatory, table_text_size=3, 
#           title_text_size=12,
#           plot_opts=list(xlab("Beta, 95% CI"), 
#                          theme(axis.title = element_text(size=12))))


### Logistic regression of IgM|IgG been positive or negative
#### IgG 
IgG.IgM_logistic_reg <- new_logistic%>%
  dplyr::filter(IgG %in% "Pos"| IgM %in% "Pos") %>%
  mutate(IgG=recode(IgG,"Pos"=1,"Neg"=0),IgM=recode(IgM,"Pos"=1,"Neg"=0))

logistic_IgG_IgM<- glm(IgG~IgG.IgM_logistic_reg$ifn_gam_quatile_cat, 
                       data = IgG.IgM_logistic_reg,family=binomial)

# IgG_reg<-IgG_logistic_reg%>%
#   mutate(IgG=ifelse(IgG=="Pos",1,0))  
# 
# IgG_reg<-IgG_logistic_reg%>%
#   mutate(IgG=recode(IgG,"Pos"=1,"Neg"=0)) 

logistic_IgG_IgM<- glm(IgG~IgG.IgM_logistic_reg$Age.category + 
                         IgG.IgM_logistic_reg$ifn_gam_quatile_cat+
                         IgG.IgM_logistic_reg$hsp_60_avrg+
                         IgG.IgM_logistic_reg$il_10_quatile_cat, 
                       data = IgG.IgM_logistic_reg,family=binomial)


### Nice Table Out Put

library(jtools) #for nice table model output
IgG_IgM.Table<-summ(logistic_IgG_IgM,
     confint = TRUE, digits = 3, 
     vifs = TRUE) # add vif to see if variance inflation factor is greater than 2
