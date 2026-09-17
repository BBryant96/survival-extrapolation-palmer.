# =============================================================================
# REFERENCES AND DATA SOURCES
# =============================================================================
#
# IPD reconstruction method:
# Guyot P, Ades AE, Ouwens MJNM, Welton NJ. Enhanced secondary analysis
# of survival data: reconstructing the data from published Kaplan-Meier
# survival curves. BMC Med Res Methodol. 2012;12:9.
# doi:10.1186/1471-2288-12-9
#
# Reconstruction software:
# Baio G. survHE: Survival Analysis for Health Economic Evaluation and
# Cost-Effectiveness Modeling. J Stat Softw. 2020;95(14):1-47.
# doi:10.18637/jss.v095.i14
#
# The survHE functions digitise() and make.ipd() were used to reconstruct
# pseudo-individual patient data from digitised Kaplan-Meier curves and
# published numbers-at-risk.
#
# ALEX Interim 1:
# Peters S, Camidge DR, Shaw AT, et al. Alectinib versus Crizotinib in
# Untreated ALK-Positive Non-Small-Cell Lung Cancer.
# N Engl J Med. 2017;377:829-838.
# doi:10.1056/NEJMoa1704795
#
# ALEX Interim 2:
# Mok T, et al. Updated overall survival and final progression-free survival
# data for patients with treatment-naive advanced ALK-positive non-small-cell
# lung cancer in the ALEX study. Ann Oncol. 2020;31(8):1056-1064.
# doi:10.1016/j.annonc.2020.04.478
#
# Digitisation software:
# Rohatgi A. WebPlotDigitizer. Available at: https://automeris.io
# Used to extract Kaplan-Meier coordinates from published ALEX trial figures.
# =============================================================================


# =============================================================================
# OS survival extrapolation — ALEX (alectinib vs crizotinib), NICE TA536
# Interim 1 data cut
# -----------------------------------------------------------------------------
# Pipeline:
#   1. Guyot IPD reconstruction from digitised KM curves + validation
# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# install.packages(c("readxl", "openxlsx", "tidyverse", "survival", "survminer",
#                    "flexsurv", "flexsurvcure", "boot", "muhaz", "broom"))
# install.packages("survHE",
#                  repos = c("https://giabaio.r-universe.dev",
#                            "https://cloud.r-project.org"))

library(survival)
library(flexsurv)        # flexsurvreg(), flexsurvspline(); also used by survHE
library(survHE)          # fit.models()
library(muhaz)           # smoothed hazard estimation
library(readxl)          # life-table import
library(broom)           # tidy() for KM step data
library(dplyr)
library(tidyr)
library(ggplot2)
library(survminer)       # ggsurvplot(), ggcoxzph()


fig_dir <- "Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)



# =============================================================================
# GUYOT IPD RECONSTRUCTION - INTERIM 1
# =============================================================================
# Workflow before this script:
#   1. Digitise each KM curve in WebPlotDigitizer; save as <arm>_OS.txt
#      (3 columns: ID, time, survival).
#   2. Record numbers at risk; save as <arm>_OS_risk.txt
#      (5 columns: Interval, Time, Lower, Upper, nrisk). "Lower"/"Upper" are the first/last digitised IDs falling in each at-risk interval.

# Reconstruct pseudo-IPD for one arm from its digitised survival + risk files.

reconstruct_ipd <- function(arm_tag, trt_label) {
  surv_file <- paste0(arm_tag, "_OS_1.txt")
  risk_file <- paste0(arm_tag, "_OS_risk_1.txt")
  km_file   <- paste0(arm_tag, "_KMdata_OS1.txt")
  ipd_file  <- paste0(arm_tag, "_IPDdata_OS1.txt")

  digitise(surv_inp   = surv_file,
           nrisk_inp  = risk_file,
           km_output  = km_file,
           ipd_output = ipd_file)

  ipd <- make.ipd(ipd_files = c(ipd_file), ctr = 1,
                  var.labs  = c("time", "event", "arm"))
  ipd$trt <- trt_label
  ipd
}

alect_ipd_os <- reconstruct_ipd("Alect", "Alectinib")
crizo_ipd_os <- reconstruct_ipd("Crizo", "Crizotinib")

save(alect_ipd_os, file = "Alect_IPD_OS1.RData")
save(crizo_ipd_os, file = "Crizo_IPD_OS1.RData")

# Bind arms; Crizotinib as reference level so the Cox HR reads alectinib vs crizotinib (matches the trial convention).

os_ipd <- bind_rows(alect_ipd_os, crizo_ipd_os)
os_ipd$trt <- factor(os_ipd$trt, levels = c("Crizotinib", "Alectinib"))


fit_km <- survfit(Surv(time, event) ~ trt, data = os_ipd)

# ---- 1.1 Reconstructed KM curves (Fig 1) ------------------------------------
km_plot <- ggsurvplot(
  fit_km,
  data              = os_ipd,
  xlab              = "Time since randomisation (months)",
  ylab              = "Overall survival probability",
  palette           = unname(trt_cols[c("Crizotinib", "Alectinib")]),
  legend.title      = "Number at risk",
  legend.labs       = c("Crizotinib", "Alectinib"),
  legend            = "bottom",
  conf.int          = TRUE,
  conf.int.alpha    = 0.15,
  censor.shape      = "|",
  censor.size       = 2.5,
  risk.table        = TRUE,
  risk.table.height = 0.22,
  risk.table.title  = "Number at risk",
  risk.table.y.text = FALSE,
  break.time.by     = 6,
  ggtheme           = theme_journal(12),
  tables.theme      = theme_journal(10) +
    theme(panel.border = element_blank(),
          panel.grid   = element_blank())
)
print(km_plot)

png(file.path(fig_dir, "Fig1_KM_Reconstructed_OS.png"),
    width = 7, height = 6, units = "in", res = 300)
print(km_plot)
dev.off()

# ---- 1.2 12-month OS --------------------------------------------------------
os_12 <- summary(fit_km, times = 12)
os_12_out <- data.frame(
  trt     = os_12$strata,
  surv_12 = os_12$surv,
  lower   = os_12$lower,
  upper   = os_12$upper
)
print(os_12_out)

# ---- 1.3 Hazard ratio -------------------------------------------------------
cox_os <- coxph(Surv(time, event) ~ trt, data = os_ipd)

get_hr <- function(cox_fit) {
  s <- summary(cox_fit)
  data.frame(
    HR      = s$coef[1, "exp(coef)"],
    lower95 = s$conf.int[1, "lower .95"],
    upper95 = s$conf.int[1, "upper .95"]
  )
}
print(get_hr(cox_os))

# ---- 1.4 Event counts (sanity check) ----------------------------------------
os_ipd %>%
  group_by(trt) %>%
  summarise(events = sum(event), n = n())


# =============================================================================
# GUYOT IPD RECONSTRUCTION - INTERIM 2
# =============================================================================

reconstruct_ipd <- function(arm_tag, trt_label) {
  digitise(surv_inp   = paste0(arm_tag, "_OS_2.txt"),
           nrisk_inp  = paste0(arm_tag, "_RT_2.txt"),
           km_output  = paste0(arm_tag, "_KMdata_OS2.txt"),
           ipd_output = paste0(arm_tag, "_IPDdata_OS2.txt"))
  
  ipd <- make.ipd(ipd_files = c(paste0(arm_tag, "_IPDdata_OS2.txt")), ctr = 1,
                  var.labs  = c("time", "event", "arm"))
  ipd$trt <- trt_label
  ipd
}

alect_ipd_os <- reconstruct_ipd("Alect", "Alectinib")
crizo_ipd_os <- reconstruct_ipd("Crizo", "Crizotinib")

save(alect_ipd_os, file = "Alect_IPD_OS2.RData")
save(crizo_ipd_os, file = "Crizo_IPD_OS2.RData")

# Bind arms; Crizotinib as reference level so the Cox HR reads alectinib vs crizotinib (matches the trial convention).

os_ipd <- bind_rows(alect_ipd_os, crizo_ipd_os)
os_ipd$trt <- factor(os_ipd$trt, levels = c("Crizotinib", "Alectinib"))
save(os_ipd, file = "IPD_OS2.RData")


# ---- 2.1 Reconstructed KM curves (Fig 2) ------------------------------------

fit_km <- survfit(Surv(time, event) ~ trt, data = os_ipd)

km_plot <- ggsurvplot(
  fit_km,
  data              = os_ipd,
  xlab              = "Time since randomisation (months)",
  ylab              = "Overall survival probability",
  palette           = unname(trt_cols[c("Crizotinib", "Alectinib")]),
  legend.title      = "Number at risk",
  legend.labs       = c("Crizotinib", "Alectinib"),
  legend            = "bottom",
  conf.int          = TRUE,
  conf.int.alpha    = 0.15,
  censor.shape      = "|",
  censor.size       = 2.5,
  risk.table        = TRUE,
  risk.table.height = 0.22,
  risk.table.title  = "Number at risk",
  risk.table.y.text = FALSE,
  break.time.by     = 6,
  ggtheme           = theme_journal(12),
  tables.theme      = theme_journal(10) +
    theme(panel.border = element_blank(),
          panel.grid   = element_blank())
)
print(km_plot)

png(file.path(fig_dir, "Fig2_KM_Reconstructed_OS.png"),
    width = 7, height = 6, units = "in", res = 300)
print(km_plot)
dev.off()


# ---- 2.2 Landmark OS (annual, 1-5 years) ------------------------------------

os_landmarks <- summary(fit_km, times = seq(12, 60, by = 12))
os_landmarks_out <- data.frame(
  trt         = os_landmarks$strata,
  time_months = os_landmarks$time,
  surv        = os_landmarks$surv,
  lower       = os_landmarks$lower,
  upper       = os_landmarks$upper
)
print(os_landmarks_out)


# ---- 2.3 Hazard ratio -------------------------------------------------------

cox_os <- coxph(Surv(time, event) ~ trt, data = os_ipd)

get_hr <- function(cox_fit) {
  s <- summary(cox_fit)
  data.frame(
    HR      = s$coef[1, "exp(coef)"],
    lower95 = s$conf.int[1, "lower .95"],
    upper95 = s$conf.int[1, "upper .95"]
  )
}
print(get_hr(cox_os))


# ---- 2.4 Event counts (sanity check) ----------------------------------------

os_ipd %>%
  group_by(trt) %>%
  summarise(events = sum(event), n = n())






