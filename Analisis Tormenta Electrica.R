## Characterization and Predictors of Mortality in Arrhythmic Storm due to Acute Ischemic Heart Disease: A Retrospective In-Hospital Cohort Analysis in a third-tier Referral Center in Mexico
## Data Analysis: Neftali Eduardo Antonio-Villa (neftalivilla@comunidad.unam.mx) 
## Latest version of Analysis 30-June-2023
## Any question regarding analysis, please contact Neftali Eduardo Antonio-Villa 

#####Library#####

library(tidyverse)
library(readxl)
library(dplyr)
library(mice)
library(haven)
library(epiR)
library(ggpubr)
library(rstatix)
library(ggthemes)
library(patchwork)
library(gtsummary)
library(data.table)
library(scales)
library(ggalluvial)
library(ggsci)
library(survival)
library(ggvenn)
library(rms)
library(survival)
library(glmnet)
library(survminer)
library(My.stepwise)
library(png)


setwd("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/Tormenta Electrica")
base <- read_excel("TE OASIS.xlsx", skip = 1)

base <- janitor::clean_names(base)

glimpse(base)
nrow(base)
base<-base%>% dplyr::mutate(id = row_number())

#####Multiple Imputation Analysis#####

base<-base%>%mutate(base, id = rownames(base))%>%mutate(id=as.numeric(id))
base2<-base%>%dplyr::select(id,bnp_pico,troponina_pico,cr_basal,ck_depi,glucosa,lactato,potasio,magnesio,fevi_total)
base2_imp<-mice::mice(base2, m=5, maxit=5,seed = 123)
base2_imp_2<-complete(base2_imp,1)
base<-base%>%dplyr::select(-c(bnp_pico,troponina_pico,cr_basal,ck_depi,glucosa,lactato,potasio,magnesio,fevi_total))%>%
  left_join(base2_imp_2,by="id")

#mice::stripplot(base2_imp, pch = 20, cex = 1.2)#ver como quedan las imputaciones junto a cada una de las variables observadas
#mice::densityplot(base2_imp)

#####Recoding Variables#####

d1<-dummies::dummy(base$cmni)
base<-cbind(base,d1)

base$PREV_CVD<-NULL
base$PREV_CVD[base$iam_1_mes==1 | base$actp_1_mes==1 | base$crvc_1==1| base$evc_1==1 | base$fa_28==1]<-1
base$PREV_CVD<-na.tools::na.replace(base$PREV_CVD,0)

base$FARMACO_CAT<-NULL
base$FARMACO_CAT[base$farmaco==0]<-1
base$FARMACO_CAT[base$farmaco>=1 & base$farmaco<=2]<-2
base$FARMACO_CAT[base$farmaco>=3]<-3
base$FARMACO_CAT<-factor(base$FARMACO_CAT,labels = c("0","1-2",">3"))

base$te_refractaria[is.na(base$te_refractaria)]<-0
base$kkiv[is.na(base$kkiv)]<-0

base$CRITICAL_STATUS<-NULL
base$CRITICAL_STATUS[base$vmi==1 | base$biac==1 | base$cap==1 |base$trfr==1 |base$mortalidad_global==1]<-1
base$CRITICAL_STATUS[is.na(base$CRITICAL_STATUS)]<-0

base$CMI_AGUDA<-NULL
base$CMI_AGUDA[base$x52==0]<-0
base$CMI_AGUDA[base$x52==1]<-1
base$CMI_AGUDA[base$x52==2]<-1

#RENAME OUTCOME
base$numero_de_vasopresores_DIC<-NULL
base$numero_de_vasopresores_DIC[base$numero_de_vasopresores>=1]<-1
base$numero_de_vasopresores_DIC[base$numero_de_vasopresores==0]<-0

#####Setting Atributes####

#Numeric
base$sodio<-as.numeric(base$sodio)
base$potasio<-as.numeric(base$potasio)
base$magnesio<-as.numeric(base$magnesio)
base$bun<-as.numeric(base$bun)
base$p_h_minimo<-as.numeric(base$p_h_minimo)
base$ast<-as.numeric(base$ast)

#Factor
base$imc<-factor(base$imc,labels = c("Underweight","Normal-Weight","Over-Weight","Obesity"))
base$cmni<-factor(base$cmni,labels = c("Non-Cardiac Comorbidities","Dilated Cardiomyopathy","Chagas Cardiomyopathy","Hypertrophic Cardiomyopathy","Other"))
base$mortalidad_global<-factor(base$mortalidad_global,labels = c("Discharged Out-of-Hospital","In-Hospital Mortality"))
base$hemodinamicamente_estable_al_ingreso[is.na(base$hemodinamicamente_estable_al_ingreso)]<-0
base$sexo<-factor(base$sexo,labels = c("Women","Men"))
base$CMI_AGUDA<-factor(base$CMI_AGUDA,labels = c("Non-Acute Coronary Syndrome","Acute Coronary Syndrome"))

#Correct Typos
base$tratamiento_icc_optimo[base$insuficiencia_cardiaca_cronica==0]<-NA
base$cardioversion[base$cardioversion==3]<-1
base$amiodarona_117[base$amiodarona_117==2]<-1
base$lidocaina[is.na(base$lidocaina)]<-0
base$beta_bloqueador_120[is.na(base$beta_bloqueador_120)]<-0
base$episodios_previos_tv_fv[is.na(base$episodios_previos_tv_fv)]<-0
base$dislipidemia[is.na(base$dislipidemia)]<-0
base$alcoholismo_actual[is.na(base$alcoholismo_actual)]<-0
base$choque_septico[is.na(base$choque_septico)]<-0
base$no_descargas_en_el_evento_de_te[base$no_descargas_en_el_evento_de_te==0]<-NA
base$patron_mixto[is.na(base$patron_mixto)]<-0
base$fv[is.na(base$fv)]<-0
base$qt_largo[base$qt_largo==1]<-1
base$taquicardia_fascicular_anterosuperior[base$taquicardia_fascicular_anterosuperior==4]<-1
base$tvp[base$tvp==0 & base$tvpns==1]<-1


base$cmni_dic<-NULL
base$cmni_dic[base$cmni=="Non-Cardiac Comorbidities"]<-0
base$cmni_dic[base$cmni!="Non-Cardiac Comorbidities"]<-1

base$mortalidad_global_rec<-NULL
base$mortalidad_global_rec[base$mortalidad_global=="In-Hospital Mortality"]<-1
base$mortalidad_global_rec[base$mortalidad_global!="In-Hospital Mortality"]<-0

###Labels##
setattr(base$sexo, "label", "Sex, (%)")
setattr(base$edad, "label", "Age, (Years)")
setattr(base$imc, "label", "BMI Categories, (%)")
setattr(base$tabaquismo_previo, "label", "Previous Smoking, (%)")
setattr(base$tabaquismo_actual, "label", "Current Smoking, (%)")
setattr(base$alcoholismo_actual, "label", "Current Alcoholism, (%)")
setattr(base$has, "label", "Hypertension, (%)")
setattr(base$dm, "label", "Diabetes, (%)")
setattr(base$dislipidemia, "label", "Dyslipidemia, (%)")
setattr(base$PREV_CVD, "label", "Previous CVD, (%)")
setattr(base$insuficiencia_cardiaca_cronica, "label", "CHF, (%)")
setattr(base$iam_1_mes, "label", "Previous Myocardial Infarction, (%)")
setattr(base$erc, "label", "CKD, (%)")
setattr(base$epoc, "label", "COPD, (%)")

table(base$cmi,base$iam_1_mes)

setattr(base$antiarritmico_previo, "label", "Use of Antiarrhythmics, (%)")
setattr(base$beta_bloqueador_32, "label", "Use of Beta-Blockers , (%)")
setattr(base$FARMACO_CAT, "label", "Number of Antiarrhythmics, (%)")


setattr(base$CMI_AGUDA, "label", "Acute Ischemic Etiology, (%)")
setattr(base$cmi, "label", "Ischemic Etiology, (%)")
setattr(base$cmni_dic, "label", "Cardiac Comorbidities, (%)")
setattr(base$cmni, "label", "Previous CVD, (%)")
setattr(base$insuficiencia_cardiaca_aguda, "label", "Acute CHF, (%)")
setattr(base$insuficiencia_cardiaca_cronica_agudizada, "label", "Chronic CHF, (%)")
setattr(base$choque_cardiogenico_67, "label", "Cardiogenic Shock, (%)")
setattr(base$edema_agudo_pulmonar, "label", "Pulmonary Edema, (%)")
setattr(base$choque_septico, "label", "Septic Shock, (%)")

setattr(base$tas, "label", "Systolic Blood Pressure, (mmHg)")
setattr(base$tad, "label", "Diastolic Blood Pressure, (mmHg)")
setattr(base$pam, "label", "Median Blood Pressure, (mmHg)")
setattr(base$fr, "label", "Respiratory-Rate, (bpm)")
setattr(base$tiempo_total_de_estancia_hospitalaria, "label", "Lenght of Hospital Stay, (Days)")
setattr(base$dia_de_estancia_en_la_que_se_presento_la_te, "label", "Time upon the onset of AS, (Days)")

setattr(base$cmni1, "label", "Dilated Cardiomyopathy, (%)")
setattr(base$cmni2, "label", "Chagas Cardiomyopathy, (%)")
setattr(base$cmni3, "label", "Hypertrophic Cardiomyopathy, (%)")
setattr(base$cmni4, "label", "Other Cardiomyopaties, (%)")
setattr(base$actp_1_mes, "label", "Percutaneous Coronary Angioplasty, (%)")
setattr(base$evc_1, "label", "Stroke, (%)")
setattr(base$fa_28, "label", "Atrial Fibrillation, (%)")

setattr(base$hemodinamicamente_estable_al_ingreso, "label", "Hemodynamically Stable, (%)")
setattr(base$fevi_total, "label", "LVEF, (%)")
setattr(base$episodios_previos_tv_fv, "label", "Previous Episodes of VT/VF, (%)")
setattr(base$killip_kimbal, "label", "Killip-Kimbal, (%)")
setattr(base$nyha, "label", "NYHA, (%)")
setattr(base$fiebre_o_infeccion_al_momento_de_la_te, "label", "Fever, (%)")

setattr(base$portador_dai, "label", "Use of ICD, (%)")
setattr(base$no_descargas_en_el_evento_de_te, "label", "Charges at ES from ICD, (n)")
setattr(base$episodio_previo_te, "label", "Previous ES, (%)")
setattr(base$te_refractaria, "label", "Refractive ES, (%)")

setattr(base$amiodarona_117, "label", "Amiodarone, (%)")
setattr(base$lidocaina, "label", "Lidocaine, (%)")
setattr(base$propafenona_119, "label", "Propafenone, (%)")
setattr(base$beta_bloqueador_120, "label", "Any Beta-Blocker, (%)")
setattr(base$tratamiento_icc_optimo, "label", "Optimal CHF Treatment, (%)")
setattr(base$cardioversion, "label", "Electrical Cardioversion, (%)")

setattr(base$bnp_pico, "label", "NT-proBNP, (pmol/l)")
setattr(base$troponina_pico, "label", "Troponin, (ng/mL)")
setattr(base$cr_basal, "label", "Serum Creatinine, (mg/dl)")
setattr(base$ck_depi, "label", "GFR by CKD-EPI, (mL/min per 1.73 m2)")
setattr(base$glucosa, "label", "Glucose, (mg/dl)")
setattr(base$lactato, "label", "Lactate, (mg/dl)")
setattr(base$potasio, "label", "Potassium, (mg/dl)")
setattr(base$magnesio, "label", "Magnesium, (mg/dl)")

setattr(base$CRITICAL_STATUS, "label", "Critical Status, (%)")
setattr(base$vmi, "label", "Invasive Mechanical Ventilation, (%)")
setattr(base$biac, "label", "Intra-Aortic Balloon Pump Therapy, (%)")
setattr(base$ecmo, "label", "Extracorporeal Membrane Oxygenation, (%)")
setattr(base$cap, "label", "Pulmonary Artery Catheter, (%)")
setattr(base$trfr, "label", "Renal Teplacement Therapy, (%)")
setattr(base$mortalidad_global, "label", "Mortality Status, (%)")

#####Descriptive Characteristics by Type of Status (Table 1)#####

base %>% 
  dplyr::select(CMI_AGUDA,
                sexo,edad,tiempo_total_de_estancia_hospitalaria,dia_de_estancia_en_la_que_se_presento_la_te,imc,
                tabaquismo_previo,alcoholismo_actual,has,dm,dislipidemia,erc,epoc,insuficiencia_cardiaca_cronica,
                PREV_CVD,iam_1_mes,actp_1_mes,crvc_1,evc_1,fa_28,cmni1,cmni2,cmni3,cmni4,
                tas,tad,pam,fr,fevi_total,killip_kimbal,nyha,fiebre_o_infeccion_al_momento_de_la_te,hemodinamicamente_estable_al_ingreso,insuficiencia_cardiaca_aguda,insuficiencia_cardiaca_cronica_agudizada,choque_cardiogenico_67,edema_agudo_pulmonar,choque_septico,
                antiarritmico_previo,FARMACO_CAT,beta_bloqueador_32,episodios_previos_tv_fv,portador_dai,no_descargas_en_el_evento_de_te,episodio_previo_te,te_refractaria,amiodarona_117,lidocaina,propafenona_119,beta_bloqueador_120,cardioversion,
                bnp_pico,troponina_pico,cr_basal,ck_depi,glucosa,lactato,
                CRITICAL_STATUS,vmi,biac,cap,trfr,mortalidad_global)%>%
  tbl_summary(by = CMI_AGUDA,missing = "no")%>%
  bold_labels()%>%
  add_overall()%>%
  add_p()%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Table_1.docx")

#####Derivation of a Risk Score for Mortality (Table 2) #####

#Baseline Model
pred.1<-coxph(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec)~ cardioversion+choque_cardiogenico_67+(bnp_pico>=10000)+te_refractaria, data=base)
summary(pred.1);BIC(pred.1)
cox.m1.ph<-cox.zph(pred.1)
survminer::ggcoxzph(cox.m1.ph)

#Diagnostic Metrics
cph(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec) ~ cardioversion+choque_cardiogenico_67+(bnp_pico>=10000)+te_refractaria, data=base)
#Points Score
round(coef(pred.1)/min(abs(coef(pred.1))),1)
round(coef(pred.1))


#Score Deployment
base$score<-base$cardioversion*2+base$choque_cardiogenico_67+as.numeric((base$bnp_pico>=10000))+base$te_refractaria
cph(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec) ~ score, data=base)

#Cox of Score Model
score.1<-coxph(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec)~ score, data=base)
summary(score.1);BIC(score.1)
cox.zph(score.1)

#Assumptions of Models
cox.m1.ph<-cox.zph(score.1)
survminer::ggcoxzph(cox.m1.ph)

#Categorical Model
base$score_cat<-NULL
base$score_cat[base$score<=1]<-1
base$score_cat[base$score>=2 & base$score<=3]<-2
base$score_cat[base$score>=4]<-3
base$score_cat<-factor(base$score_cat)
table(base$score_cat,base$mortalidad_global_rec)

#Cox of Score Model (Categorical)
score.2<-coxph(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec)~ score_cat, data=base)
summary(score.2);BIC(score.2)
cox.m3.ph<-cox.zph(score.2)
survminer::ggcoxzph(cox.m3.ph)


cph(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec) ~ score_cat, data=base)

#####Types of EKG Patherns (Figure 3)######

Figure3A.df<-base%>%
  dplyr::select(tvms,fv,tvp)%>%
  mutate("tvms"=as.logical(tvms),
         "fv"=as.logical(fv),
         "tvp"=as.logical(tvp))

names(Figure3A.df)<-c("Sustained Monomorphic VT","Ventricular Fibrilation","Polymorphic VT")

Figure3A<-ggvenn(Figure3A.df,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 4)+labs(title = "EKG Patterns in Arrhythmic Storm")

Figure3B<-base %>% 
  dplyr::select(105:113,qt_largo,taquicardia_fascicular_anterosuperior)%>%
  summarise(across(everything(), list(sum),na.rm=T))%>%
  tidyr::gather(group, value)%>%
  mutate(perc.num=value / nrow(base)*100,
         perc.lab = paste(paste0("(",round(value / nrow(base)*100,2),"%",")")))%>%
  ggplot(aes(reorder(group, -perc.num), perc.num,fill=group)) +
  geom_col(fill=palette_pander(11),col="black",alpha=0.75)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()+
  scale_x_discrete(labels= c("Right Bundle Branch Block",
                             "Left Bundle Branch Block",
                             "Complete Heart Block",
                             "Atrial Fibrilation",
                             "First Degree Heart Block",
                             "Left Anterior Fascicular Block",
                             "Idiopathic Fascicular LV Tachycardia",
                             "Second Degree Heart Block",
                             "Brugada syndrome",
                             "Long QT syndrome",
                             "Left Posterior Fascicular Block"))+
  xlab("")+
  ylab("Percentage, (%)")+
  theme_classic()+
  labs(fill="",
       title ="Arrhythmias Before Arrhythmic Storm",
       subtitle = "n=101")+
  geom_text(aes(label=perc.lab), vjust = 0.5,hjust = -0.2, color="black", size=4)+
  scale_y_continuous(limits = c(0, 40))+
  guides(fill="none")


Figure3<-ggarrange(Figure3B,Figure3A,ncol = 2,nrow = 1,labels = LETTERS[1:2],widths=c(1,0.75))

ggsave(file = "Figure3.pdf", 
       Figure3,
       bg = "transparent",
       width = 35, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)


#####Derivation of a Risk Score for Mortality (Figure 4)#####

#Overall Population (Stepwise Model)
stepwiseCox(formula=Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec) ~ 
              crvc_1+erc+beta_bloqueador_32+antiarritmico_previo+has+choque_cardiogenico_67+choque_septico+fr+edad+
              has+hemodinamicamente_estable_al_ingreso+episodios_previos_tv_fv+portador_dai+episodio_previo_te+te_refractaria+brdhh+bloqueo_de_fasciculo_anterior+cardioversion+bnp_pico+ck_depi+CMI_AGUDA,
            data=base,
            selection="bidirection",
            select="SBC",
            sls=0.05,
            method="breslow")

#Overall Population (Final Model)
m1<-coxph(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec)~ cardioversion+choque_cardiogenico_67+(bnp_pico>=10000)+te_refractaria, 
          data=base,method = "breslow")
summary(m1);BIC(m1);car::vif(m1)
cox.zph(m1)

#Plot Model
Figure4A<-jtools::plot_summs(m1,
                             colors = c(ggsci::pal_aaas()(1)),
                             model.names = c("All-Sample"),
                             exp=T,coefs = c("Electrical Cardioversion, (%)"="cardioversion",
                                             "Cardiogenic Shock, (%)"="choque_cardiogenico_67",
                                             "NT-proBNP >10,000 pg/ml, (%)"="bnp_pico >= 10000TRUE",
                                             "Refractive AS, (%)"="te_refractaria"))+
  xlab("Adjusted Cox Regression Model: \nHazard Ratio (HR, 95% CI)")+
  ylab("")+
  ggtitle("")+
  theme_classic() +
  theme(axis.line=element_blank(),panel.grid.minor = element_blank(), axis.ticks.y = element_blank())+
  guides(colour = guide_legend(nrow = 3))+
  scale_x_log10()+
  labs(title ="Predictors of In-Hospital Mortality")+
  geom_text(inherit.aes = FALSE, aes(x = 3.117, y = 1.1, label = "HR: 2.117 (95% CI: 1.188-3.771, p=0.011)", vjust = -0.5))+
  geom_text(inherit.aes = FALSE, aes(x = 3.330, y = 2.1, label = "HR: 2.330 (95% CI: 1.277-4.251, p=0.005)", vjust = -0.5))+
  geom_text(inherit.aes = FALSE, aes(x = 3.353, y = 3.1, label = "HR: 2.353 (95% CI: 1.014-5.463, p=0.046)", vjust = -0.5))+
  geom_text(inherit.aes = FALSE, aes(x = 6.341, y = 4.1, label = "HR: 6.341 (95% CI: 2.632-15.275, p<0.001)", vjust = -0.5))

#Table of Points
img1 <- readPNG("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/Tormenta Electrica/Figure_4A.png")

Figure4B<-ggplot() + 
  background_image(img1) +
  theme(plot.margin = margin(t=0.75, l=0.75, r=0.75, b=0.75, unit = "cm"))+
  ggtitle("Mortality Risk Score")

#Predicted Probability
base$risk_probability<-(1-predict(score.1,type="survival"))*100
base$survival_probability<-(predict(score.1,type="survival"))*100

#Evaluation of Risk Probability by Categories
base$score_cat<-factor(base$score_cat,labels = c("Mid-Risk","High-Risk","Very\nHigh-Risk"))
my_comparisons <- list( c("Mid-Risk", "High-Risk"), c("Mid-Risk", "Very\nHigh-Risk"),c("High-Risk", "Very\nHigh-Risk"))

Figure4C<-ggplot(base,aes(score_cat,survival_probability,fill=score_cat))+
  geom_boxplot()+
  theme_classic()+
  xlab("")+
  ylab("Probability, (%)")+
  stat_compare_means(comparisons = my_comparisons)+
  scale_fill_manual(values = c("#2dc937", "#db7b2b","#cc3232"))+
  labs(fill="Score\nClasification")+
  ggtitle("Predicted Survival Probability")+
  theme(legend.position = "top")


#Kaplan Meier Curve Analisis
mod1_km<-survfit(Surv(dias_de_estancia_hospitalaria_a_partir_de_la_te, mortalidad_global_rec) ~factor(score_cat), data = base)
KM_fig1<-ggsurvplot(mod1_km, data = base, size = 1,palette =c("#2dc937", "#db7b2b","#cc3232"),conf.int = T,
                    risk.table = T,
                    ggtheme = theme_classic(),
                    xlab="Time since the onset of AS, (Days)",
                    ylab="Survival Probability, (%)",
                    title="Kaplan-Meier Survival Probability",
                    legend.labs = c("Mid-Risk",
                                    "High-Risk",
                                    "Very High-Risk"),
                    xlim = c(0,28),
                    ylim= c(0,1.0),
                    break.y.by= c(0.1),
                    break.x.by= c(7),
                    pval = TRUE, 
                    pval.method = TRUE,
                    log.rank.weights = "1", 
                    pval.method.size = 3,
                    pval.coord = c(2, 0.10),
                    pval.method.coord = c(2, 0.05))+
  theme_survminer(base_size = 9,base_family = "Arial")+
  guides(colour = guide_legend(nrow = 1))

KM_fig1 <-KM_fig1 + theme_survminer(base_size = 10,
                                    base_family = "Helvetica",
                                    font.x = c(10, "plain" ), 
                                    font.y = c(10, "plain"),
                                    font.main = c(10, "plain"),
                                    font.caption = c(8, "plain"), 
                                    font.legend = c(8, "plain"),
                                    font.tickslab = c(8, "plain"))

Figure4D<-ggarrange(KM_fig1$plot, KM_fig1$table, heights = c(2, 0.7),
               ncol = 1, nrow = 2);Figure4D

Figure4<-ggarrange(Figure4A,Figure4B,Figure4C,Figure4D,ncol = 2,nrow = 2,labels = c(LETTERS[1:4]))

ggsave(Figure4, 
       file="Figure4.pdf", 
       bg="transparent",
       width=32.5, height=25, 
       units=c("cm"), dpi=600, limitsize = FALSE)


#####Pharmacological and Laboratories stratified by ACS (Supplementary Table 2)#####

base %>% 
  dplyr::select(CMI_AGUDA,31:50,175:192)%>%
  dplyr::select(-c(farmaco,tratamiento_icc_optimo,id))%>%
  tbl_summary(by = CMI_AGUDA,
              missing = "no")%>%
  bold_labels()%>%
  add_overall()%>%
  add_p()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Supplementary_Table_2.docx")

#####Descriptive Characteristics by Mortality Status (Supplementary Table 3)#####

base %>% 
  dplyr::select(sexo,edad,tiempo_total_de_estancia_hospitalaria,dia_de_estancia_en_la_que_se_presento_la_te,imc,
    tabaquismo_previo,alcoholismo_actual,has,dm,dislipidemia,erc,epoc,insuficiencia_cardiaca_cronica,
    PREV_CVD,iam_1_mes,actp_1_mes,crvc_1,evc_1,fa_28,cmni1,cmni2,cmni3,cmni4,
    tas,tad,pam,fr,fevi_total,killip_kimbal,nyha,fiebre_o_infeccion_al_momento_de_la_te,hemodinamicamente_estable_al_ingreso,insuficiencia_cardiaca_aguda,insuficiencia_cardiaca_cronica_agudizada,choque_cardiogenico_67,edema_agudo_pulmonar,choque_septico,
    antiarritmico_previo,FARMACO_CAT,beta_bloqueador_32,episodios_previos_tv_fv,portador_dai,no_descargas_en_el_evento_de_te,episodio_previo_te,te_refractaria,amiodarona_117,lidocaina,propafenona_119,beta_bloqueador_120,cardioversion,
    bnp_pico,troponina_pico,cr_basal,ck_depi,glucosa,lactato,
    CRITICAL_STATUS,vmi,biac,cap,trfr,mortalidad_global,CMI_AGUDA)%>%
  tbl_summary(by = mortalidad_global,missing = "no")%>%
  bold_labels()%>%
  add_overall()%>%
  add_p()%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Supplementary_Table_3.docx")


#####Univariate Associatied Variables by Mortality Status (Supplementary Table 4)#####

base$killip_kimbal<-factor(base$killip_kimbal)
base$nyha<-factor(base$nyha)

setattr(base$brdhh, "label", "Right Bundle Branch Block, (%)")
setattr(base$brihh, "label", "Left Bundle Branch Block, (%)")
setattr(base$bav_completo, "label", "Complete Heart Block, (%)")
setattr(base$fa_114, "label", "Atrial Fibrilation, (%)")
setattr(base$bav_1er_grado, "label", "First Degree Heart Block, (%)")
setattr(base$bav_2do, "label", "Second Degree Heart Block, (%)")
setattr(base$bloqueo_de_fasciculo_anterior, "label", "Left Anterior Fascicular Block, (%)")
setattr(base$bloqueo_de_fasciculo_psoterior, "label", "Left Posterior Fascicular Block, (%)")
setattr(base$brugada, "label", "Brugada syndrome, (%)")
setattr(base$qt_largo, "label", "Long QT syndrome, (%)")
setattr(base$taquicardia_fascicular_anterosuperior, "label", "Idiopathic Fascicular LV Tachycardia, (%)")
setattr(base$tvms, "label", "Sustained Monomorphic VT, (%)")
setattr(base$fv, "label", "Ventricular Fibrilation, (%)")
setattr(base$tvp, "label", "Polymorphic VT, (%)")

#Overall Population
sup.tab.4.1<-tbl_uvregression(
  base,
  method=coxph,
  y = Surv(time = dias_de_estancia_hospitalaria_a_partir_de_la_te, event = base$mortalidad_global=="In-Hospital Mortality"),
  exponentiate = TRUE,
  hide_n = TRUE,
  include = c(sexo,edad,tiempo_total_de_estancia_hospitalaria,dia_de_estancia_en_la_que_se_presento_la_te,imc,
              tabaquismo_previo,alcoholismo_actual,has,dm,dislipidemia,erc,insuficiencia_cardiaca_cronica,
              PREV_CVD,iam_1_mes,actp_1_mes,crvc_1,evc_1,fa_28,cmni1,cmni2,cmni3,cmni4,
              tas,tad,pam,fr,fevi_total,killip_kimbal,nyha,fiebre_o_infeccion_al_momento_de_la_te,hemodinamicamente_estable_al_ingreso,insuficiencia_cardiaca_aguda,insuficiencia_cardiaca_cronica_agudizada,choque_cardiogenico_67,edema_agudo_pulmonar,choque_septico,
              antiarritmico_previo,FARMACO_CAT,beta_bloqueador_32,episodios_previos_tv_fv,portador_dai,no_descargas_en_el_evento_de_te,episodio_previo_te,te_refractaria,amiodarona_117,lidocaina,beta_bloqueador_120,cardioversion,
              bnp_pico,troponina_pico,cr_basal,ck_depi,glucosa,lactato,mortalidad_global,CMI_AGUDA,
              tvms,fv,tvp,105:113,qt_largo,taquicardia_fascicular_anterosuperior))%>%
  bold_labels()%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label) )

#Patients with ACS

sup.tab.4.2<-tbl_uvregression(
  base%>%dplyr::filter(CMI_AGUDA=="Acute Coronary Syndrome"),
  method=coxph,
  y = Surv(time = dias_de_estancia_hospitalaria_a_partir_de_la_te, event = mortalidad_global=="In-Hospital Mortality"),
  exponentiate = TRUE,
  hide_n = TRUE,
  include = c(sexo,edad,tiempo_total_de_estancia_hospitalaria,dia_de_estancia_en_la_que_se_presento_la_te,imc,
              tabaquismo_previo,alcoholismo_actual,has,dm,dislipidemia,erc,insuficiencia_cardiaca_cronica,
              PREV_CVD,iam_1_mes,actp_1_mes,crvc_1,evc_1,fa_28,cmni1,cmni2,cmni3,cmni4,
              tas,tad,pam,fr,fevi_total,killip_kimbal,nyha,fiebre_o_infeccion_al_momento_de_la_te,hemodinamicamente_estable_al_ingreso,insuficiencia_cardiaca_aguda,insuficiencia_cardiaca_cronica_agudizada,choque_cardiogenico_67,edema_agudo_pulmonar,choque_septico,
              antiarritmico_previo,FARMACO_CAT,beta_bloqueador_32,episodios_previos_tv_fv,portador_dai,no_descargas_en_el_evento_de_te,episodio_previo_te,te_refractaria,amiodarona_117,lidocaina,beta_bloqueador_120,cardioversion,
              bnp_pico,troponina_pico,cr_basal,ck_depi,glucosa,lactato,mortalidad_global,CMI_AGUDA,
              tvms,fv,tvp,105:113,qt_largo,taquicardia_fascicular_anterosuperior))%>%
  bold_labels()%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))

#Patients without ACS

sup.tab.4.3<-tbl_uvregression(
  base%>%dplyr::filter(CMI_AGUDA!="Acute Coronary Syndrome"),
  method=coxph,
  y = Surv(time = dias_de_estancia_hospitalaria_a_partir_de_la_te, event = mortalidad_global=="In-Hospital Mortality"),
  exponentiate = TRUE,
  hide_n = TRUE,
  include = c(sexo,edad,tiempo_total_de_estancia_hospitalaria,dia_de_estancia_en_la_que_se_presento_la_te,imc,
              tabaquismo_previo,alcoholismo_actual,has,dm,dislipidemia,erc,insuficiencia_cardiaca_cronica,
              PREV_CVD,iam_1_mes,actp_1_mes,crvc_1,evc_1,fa_28,cmni1,cmni2,cmni3,cmni4,
              tas,tad,pam,fr,fevi_total,killip_kimbal,nyha,fiebre_o_infeccion_al_momento_de_la_te,hemodinamicamente_estable_al_ingreso,insuficiencia_cardiaca_aguda,insuficiencia_cardiaca_cronica_agudizada,choque_cardiogenico_67,edema_agudo_pulmonar,choque_septico,
              antiarritmico_previo,FARMACO_CAT,beta_bloqueador_32,episodios_previos_tv_fv,portador_dai,no_descargas_en_el_evento_de_te,episodio_previo_te,te_refractaria,amiodarona_117,lidocaina,beta_bloqueador_120,cardioversion,
              bnp_pico,troponina_pico,cr_basal,ck_depi,glucosa,lactato,mortalidad_global,CMI_AGUDA,
              tvms,fv,tvp,105:113,qt_largo,taquicardia_fascicular_anterosuperior))%>%
  bold_labels()%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))

tbl_merge(
  tbls = list(sup.tab.4.1,sup.tab.4.2,sup.tab.4.3),
  tab_spanner = c("**Overall Population**", "**Acute Coronary Syndrome**","**Non-ACS**"))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Supplementary_Table_4.docx")


#####Types of EKG Patherns Stratified by Acute Myocardial Infarction (Suplementary Figure 1)######

Sup.Figure3A.1.df<-base%>%
  dplyr::filter(CMI_AGUDA==0)%>%
  dplyr::select(tvms,fv,tvp)%>%
  mutate("tvms"=as.logical(tvms),
         "fv"=as.logical(fv),
         "tvp"=as.logical(tvp))

names(Sup.Figure3A.1.df)<-c("Sustained Monomorphic VT","Ventricular Fibrilation","Polymorphic VT")

Sup.Figure3A.1<-ggvenn(Sup.Figure3A.1.df,
                       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
                       stroke_size = 0.5, set_name_size = 4)+labs(title = "EKG Patterns of AS\n(Non-ACS)")

Sup.Figure3A.2.df<-base%>%
  dplyr::filter(CMI_AGUDA==1)%>%
  dplyr::select(tvms,fv,tvp)%>%
  mutate("tvms"=as.logical(tvms),
         "fv"=as.logical(fv),
         "tvp"=as.logical(tvp))

names(Sup.Figure3A.2.df)<-c("Sustained Monomorphic VT","Ventricular Fibrilation","Polymorphic VT")

Sup.Figure3A.2<-ggvenn(Sup.Figure3A.2.df,
                       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
                       stroke_size = 0.5, set_name_size = 4)+labs(title = "EKG Patterns of AS\n(Acute Coronary Syndrome)")

Sup.Figure3A<-ggarrange(Sup.Figure3A.1,Sup.Figure3A.2,ncol = 2,nrow = 1)

#Stratification by Type of Coronary Syndrome

Sup.Figure3B<-base %>% 
  dplyr::select(CMI_AGUDA,105:113,qt_largo,taquicardia_fascicular_anterosuperior)%>%
  dplyr::group_by(CMI_AGUDA)%>%
  summarise(across(everything(), list(sum),na.rm=T))%>%
  gather(group, value,-CMI_AGUDA)%>%
  mutate(n=rep(c(55,46),11),
         perc.num=value / n*100,
         perc.lab = paste(paste0("(",round(value / n*100,2),"%",")")))%>%
  ggplot(aes(reorder(group, -perc.num), perc.num,fill=factor(CMI_AGUDA),group=CMI_AGUDA)) +
  geom_col(col="black",alpha=0.75,position = "dodge")+
  scale_fill_nejm()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()+
  scale_x_discrete(labels= c("Right Bundle Branch Block",
                             "Left Bundle Branch Block",
                             "Complete Heart Block",
                             "Atrial Fibrilation",
                             "First Degree Heart Block",
                             "Left Anterior Fascicular Block",
                             "Second Degree Heart Block",
                             "Brugada syndrome",
                             "Long QT syndrome",
                             "Anterior-Superior Fascicular Tachycardia",
                             "Left Posterior Fascicular Block"))+
  xlab("")+
  ylab("Percentage, (%)")+
  theme_classic()+
  labs(fill="",
       title ="ECG Patterns Before Arrhythmic Storm",
       subtitle = "n=101")+
  geom_text(aes(label=perc.lab), position = position_dodge(width = 0.9), vjust = 0.5,hjust = -0.2, color="black", size=4)+
  scale_y_continuous(limits = c(0, 40))+
  theme(legend.position = "top")

chisq.test(table(base$brdhh,base$CMI_AGUDA))
chisq.test(table(base$brihh,base$CMI_AGUDA))
chisq.test(table(base$bav_completo,base$CMI_AGUDA))
chisq.test(table(base$bloqueo_de_fasciculo_anterior,base$CMI_AGUDA))
chisq.test(table(base$bloqueo_de_fasciculo_psoterior,base$CMI_AGUDA))
chisq.test(table(base$bav_2do,base$CMI_AGUDA))
chisq.test(table(base$bav_1er_grado,base$CMI_AGUDA))
chisq.test(table(base$brugada,base$CMI_AGUDA))
chisq.test(table(base$fa_114,base$CMI_AGUDA))
chisq.test(table(base$tvms,base$CMI_AGUDA))
chisq.test(table(base$fv,base$CMI_AGUDA))
chisq.test(table(base$tvp,base$CMI_AGUDA))

Sup.Figure3<-ggarrange(Sup.Figure3B,Sup.Figure3A,ncol = 2,nrow = 1,labels = c("A","B"),widths=c(1,1))

ggsave(file = "Sup.Figure1.pdf", 
       Sup.Figure3,
       bg = "transparent",
       width = 40, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)



#####Assumtion Models of Mortality Score (Supplementary Figure 2)#####
#####Causes of Death and Predictors or Mortality (Supplementary Figure 3)#####

#Causes of Death

Sup.Figure.3A<-base %>% 
  dplyr::filter(mortalidad_global=="In-Hospital Mortality")%>%
  dplyr::select(choque_cardiogenico_142,eap,iam,complicacion_mecanica_post_iam,
                sepsis,evc,falla_renal,choque_hipovolemico,tep_154,neumonia)%>%
  summarise(across(everything(), list(sum),na.rm=T))%>%
  tidyr::gather(group, value)%>%
  mutate(perc.num=value / nrow(base)*100,
         perc.lab = paste(paste0("(",round(value / nrow(base)*100,2),"%",")")))%>%
  ggplot(aes(reorder(group, -perc.num), perc.num,fill=group)) +
  geom_col(fill=palette_pander(10),col="black",alpha=0.75)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()+
  scale_x_discrete(labels= c("Cardiogenic Shock",
                             "Death by Myocardial Infarction",
                             "Acute Kidney Faillure",
                             "Sepsis",
                             "Pneumonia",
                             "Mechanical Complication by Myocardial Infarction",
                             "Hypovolemic Shock",
                             "Acute Pulmonary Edema",
                             "Stroke",
                             "Pulmonary Embolism"))+
  xlab("")+
  ylab("Percentage, (%)")+
  theme_classic()+
  labs(fill="",
       title ="Causes of In-Hospital Mortality",
       subtitle = "n=49")+
  geom_text(aes(label=perc.lab), vjust = 0.5,hjust = -0.2, color="black", size=4)+
  scale_y_continuous(limits = c(0, 100))+
  guides(fill="none")


Sup.Figure.3B<-base %>% 
  dplyr::filter(mortalidad_global=="In-Hospital Mortality")%>%
  dplyr::select(CMI_AGUDA,choque_cardiogenico_142,eap,iam,complicacion_mecanica_post_iam,
                sepsis,evc,falla_renal,choque_hipovolemico,tep_154,neumonia)%>%
  dplyr::group_by(CMI_AGUDA)%>%
  summarise(across(everything(), list(sum),na.rm=T))%>%
  tidyr::gather(group, value,-CMI_AGUDA)%>%
  mutate(n=rep(c(21,28),10),
         perc.num=value / n*100,
         perc.lab = paste(paste0("(",round(value / n*100,2),"%",")")))%>%
  ggplot(aes(reorder(group, -perc.num), perc.num,fill=factor(CMI_AGUDA),group=CMI_AGUDA)) +
  geom_col(col="black",alpha=0.75,position = "dodge")+
  scale_fill_nejm()+
  coord_flip()+
  scale_x_discrete(labels= c("Cardiogenic Shock",
                             "Death by Myocardial Infarction",
                             "Acute Kidney Faillure",
                             "Sepsis",
                             "Pneumonia",
                             "Mechanical Complication by Myocardial Infarction",
                             "Acute Pulmonary Edema",
                             "Pulmonary Embolism",
                             "Hypovolemic Shock",
                             "Stroke"))+
  xlab("")+
  ylab("Percentage, (%)")+
  theme_classic()+
  labs(fill="",
       title ="Causes of In-Hospital Mortality",
       subtitle = "n=49")+
  geom_text(aes(label=perc.lab), position = position_dodge(width = 0.9), vjust = 0.5,hjust = -0.2, color="black", size=4)+
  scale_y_continuous(limits = c(0, 100))+
  theme(legend.position = "top")

Sup.Figure.3<-ggarrange(Sup.Figure.3A,Sup.Figure.3B,ncol = 2,nrow = 1,labels = LETTERS[1:2],widths = c(0.75,1.0))

ggsave(file = "Sup.Figure3.pdf", 
       Sup.Figure.3,
       bg = "transparent",
       width = 40, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)

#####Mortality Rate#####

#Overall population
sm0<-base %>% filter(dias_de_estancia_hospitalaria_a_partir_de_la_te>=0) %>% group_by(mortalidad_global=="In-Hospital Mortality") %>%
  summarise(cases=n(), time=sum(dias_de_estancia_hospitalaria_a_partir_de_la_te)/7)%>%na.omit()

time<-sum(sm0$time)
cases<-sm0$cases[sm0$`mortalidad_global == "In-Hospital Mortality"`==TRUE]
ncas <- cases; ntar <- time
tmp <- as.matrix(cbind(ncas, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 1000, design = 1, 
         conf.level = 0.95) * 100

#With ACS

sm0<-base %>% filter(dias_de_estancia_hospitalaria_a_partir_de_la_te>=0) %>% 
  filter(CMI_AGUDA=="Acute Coronary Syndrome") %>% 
  group_by(mortalidad_global=="In-Hospital Mortality") %>%
  summarise(cases=n(), time=sum(dias_de_estancia_hospitalaria_a_partir_de_la_te)/7)%>%na.omit()

time<-sum(sm0$time)
cases<-sm0$cases[sm0$`mortalidad_global == "In-Hospital Mortality"`==TRUE]
ncas <- cases; ntar <- time
tmp <- as.matrix(cbind(ncas, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 1000, design = 1, 
         conf.level = 0.95) * 100

#Without ACS

sm0<-base %>% filter(dias_de_estancia_hospitalaria_a_partir_de_la_te>=0) %>% 
  filter(CMI_AGUDA!="Acute Coronary Syndrome") %>% 
  group_by(mortalidad_global=="In-Hospital Mortality") %>%
  summarise(cases=n(), time=sum(dias_de_estancia_hospitalaria_a_partir_de_la_te)/7)%>%na.omit()

time<-sum(sm0$time)
cases<-sm0$cases[sm0$`mortalidad_global == "In-Hospital Mortality"`==TRUE]
ncas <- cases; ntar <- time
tmp <- as.matrix(cbind(ncas, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 1000, design = 1, 
         conf.level = 0.95) * 100
