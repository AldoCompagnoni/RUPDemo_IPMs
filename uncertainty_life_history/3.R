

# Life history pipeline -------------------------------------------------------
source("uncertainty_life_history/life_history_pipeline.R")

# Models ----------------------------------------------------------------------
models <- c(
  "adler_2007_ks/data/ange/adler07_ange_mean_ipmr.rds",
  "adler_2007_ks/data/arlo/adler07_arlo_mean_ipmr.rds",
  "adler_2007_ks/data/bocu/adler07_bocu_mean_ipmr.rds",
  "adler_2007_ks/data/bogr/adler07_bogr_mean_ipmr.rds",
  "adler_2007_ks/data/bohi/adler07_bohi_mean_ipmr.rds",
  "adler_2007_ks/data/buda/adler07_buda_mean_ipmr.rds",
  "adler_2007_ks/data/pavi/adler07_pavi_mean_ipmr.rds",
  "adler_2007_ks/data/scsc/adler07_scsc_mean_ipmr.rds",
  "adler_2007_ks/data/spas/adler07_spas_mean_ipmr.rds",
  "adler_2007_ks/data/spcr/adler07_spcr_mean_ipmr.rds",
  "anderson_2016_az/data/boer/anderson16az_boer_mean_ipmr.rds",
  "anderson_2016_az/data/boro/anderson16az_boro_mean_ipmr.rds",
  "anderson_2016_az/data/dica/anderson16az_dica_mean_ipmr.rds",
  "anderson_2016_az/data/heco/anderson16az_heco_mean_ipmr.rds",
  "anderson_2016_az/data/hibe/anderson16az_hibe_mean_ipmr.rds",
  "anderson_2016_az/data/lyph/anderson16az_lyph_mean_ipmr.rds",
  "anderson_2016_mt/data/arfr/anderson16mt_arfr_mean_ipmr.rds",
  "anderson_2016_mt/data/boda/anderson16mt_boda_mean_ipmr.rds",
  "anderson_2016_mt/data/bogr/anderson16mt_bogr_mean_ipmr.rds",
  "anderson_2016_mt/data/cafi/anderson16mt_cafi_mean_ipmr.rds",
  "anderson_2016_mt/data/heco/anderson16mt_heco_mean_ipmr.rds",
  "anderson_2016_mt/data/pose/anderson16mt_pose_mean_ipmr.rds",
  "anderson_2016_mt/data/scpa/anderson16mt_scpa_mean_ipmr.rds",
  #"anderson_2016_mt/data/spco/anderson16mt_spco_mean_ipmr.rds",
  "christensen_2021_nm/data/arad/christensen21_arad_mean_ipmr.rds",
  "christensen_2021_nm/data/boer/christensen21_boer_mean_ipmr.rds",
  "christensen_2021_nm/data/dapu/christensen21_dapu_mean_ipmr.rds",
  "christensen_2021_nm/data/muar/christensen21_muar_mean_ipmr.rds",
  "christensen_2021_nm/data/paob/christensen21_paob_mean_ipmr.rds",
  "christensen_2021_nm/data/plmu/christensen21_plmu_mean_ipmr.rds",
  "christensen_2021_nm/data/prgl/christensen21_prgl_mean_ipmr.rds",
  "christensen_2021_nm/data/scbr/christensen21_scbr_mean_ipmr.rds",
  "christensen_2021_nm/data/spai/christensen21_spai_mean_ipmr.rds",
  "christensen_2021_nm/data/spfl/christensen21_spfl_mean_ipmr.rds",
  "chu_2013_co/data/arlo/chu13_arlo_mean_ipmr.rds",
  "chu_2013_co/data/buda/chu13_buda_mean_ipmr.rds",
  "chu_2013_co/data/scpa/chu13_scpa_mean_ipmr.rds",
  "chu_2013_co/data/sihy/chu13_sihy_mean_ipmr.rds",
  "chu_2013_co/data/spcr/chu13_spcr_mean_ipmr.rds",
  "moore_2021_az/data/anpa/moore21_anpa_mean_ipmr.rds",
  "moore_2021_az/data/anro/moore21_anro_mean_ipmr.rds",
  "moore_2021_az/data/bogr/moore21_bogr_mean_ipmr.rds",
  "moore_2021_az/data/cage/moore21_cage_mean_ipmr.rds",
  "moore_2021_az/data/elel/moore21_elel_mean_ipmr.rds",
  "moore_2021_az/data/eltr/moore21_eltr_mean_ipmr.rds",
  "moore_2021_az/data/erfe/moore21_erfe_mean_ipmr.rds",
  "moore_2021_az/data/fear/moore21_fear_mean_ipmr.rds",
  "moore_2021_az/data/koma/moore21_koma_mean_ipmr.rds",
  "moore_2021_az/data/mumo/moore21_mumo_mean_ipmr.rds",
  "moore_2021_az/data/muri/moore21_muri_mean_ipmr.rds",
  "moore_2021_az/data/mutr/moore21_mutr_mean_ipmr.rds",
  "moore_2021_az/data/muwr/moore21_muwr_mean_ipmr.rds",
  "moore_2021_az/data/poco/moore21_poco_mean_ipmr.rds",
  "moore_2021_az/data/pofe/moore21_pofe_mean_ipmr.rds",
  "moore_2021_az/data/popr/moore21_popr_mean_ipmr.rds",
  "moore_2021_az/data/spin/moore21_spin_mean_ipmr.rds",
  "zachmann_2016_id/data/artr/zachmann16_artr_mean_ipmr.rds",
  "zachmann_2016_id/data/heco/zachmann16_heco_mean_ipmr.rds",
  "zachmann_2016_id/data/pose/zachmann16_pose_mean_ipmr.rds",
  "zachmann_2016_id/data/pssp/zachmann16_pssp_mean_ipmr.rds")

life_history_df <- run_lh_pipeline(models, overwrite = FALSE)

life_history_df





  