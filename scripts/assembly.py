#%%
import os
try:
    os.chdir("/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/igo/assay/assembly")
    print(os.getcwd())
except:
	pass

#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotnine
import pandas as pd
from glob import glob
from pprint import pprint
from plotnine import *

#%%
save_dir = "/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/igo/assay/assembly/results"
staff_oral  = "/hwfssz5/ST_META/P18Z10200N0127_MA/zhujie/sgb_oral"
yunnan_oral = "/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/yunnan_oral"
pub1_oral   = "/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/pub_oral"
pub2_oral   = "/hwfssz1/ST_META/F16ZQSB1SY2794_RO_GRJ/RO2"

#%%
def sample_site(row):
    if row["sample_id"].endswith("saliva"):
        return "saliva"
    if row["sample_id"].endswith("togue"):
        return "tongue"
    if row["sample_id"].endswith("tongue"):
        return "tongue"
    else:
        return "other"

def sample_site2(row):
    if row["bin_id"].split(".")[0].endswith("saliva"):
        return "saliva"
    if row["bin_id"].split(".")[0].endswith("togue"):
        return "tongue"
    if row["bin_id"].split(".")[0].endswith("tongue"):
        return "tongue"
    else:
        return "other"

def get_status_df(asm_tsv, project, is_bin=False):
    choose_items = ["sample_id",
                    "n_scaffolds",
                    "scaf_bp", 
                    "scaf_N50",
                    "scaf_L50",
                    "scaf_N90",
                    "scaf_L90", 
                    "scaf_max",
                    "scaf_n_gt50K",
                    "scaf_pct_gt50K",
                    "gc_avg",
                    "gc_std",
                    "project"]
    df = pd.read_csv(asm_tsv, sep='\t')\
           .rename(columns={"scaf_N50": "scaf_L50",
                            "scaf_L50": "scaf_N50",
                            "scaf_N90": "scaf_L90",
                            "scaf_L90": "scaf_N90"})\
           .assign(project=project)
    df = df.assign(sample_id=df.filename.apply(lambda x: x.split("/")[-1].split(".")[0]))
    if is_bin:
        df = df.assign(bin_id=df.filename.apply(lambda x: x.split("/")[-1].rpartition(".")[0]))
        df = df.assign(bin_id_2=df.bin_id.apply(lambda x: x.split(".")[0] + "." + x.split(".")[-1]))
        return df.loc[:, ["bin_id", "bin_id_2"] + choose_items]
    else:
        return df.loc[:, choose_items]

#%%
#### staff oral
staff_df_m0 = []
for i in glob(staff_oral + "/assay/03.assembly/status/m0/*.tsv"):
    staff_df_m0 = staff_df_m0 + [get_status_df(i, "staff_oral")]
for i in glob(staff_oral + "/assay/03.assembly/status/m0_1/*.tsv"):
    staff_df_m0 = staff_df_m0 + [get_status_df(i, "staff_oral")]
for i in glob(staff_oral + "/assay/03.assembly/status/m0_2/*.tsv"):
    staff_df_m0 = staff_df_m0 + [get_status_df(i, "staff_oral")]    
staff_asm_df_m0 = pd.concat(staff_df_m0)

staff_df_m1000 = []
for i in glob(staff_oral + "/assay/03.assembly/status/m1000/*.tsv"):
    staff_df_m1000 = staff_df_m1000 + [get_status_df(i, "staff_oral")]
for i in glob(staff_oral + "/assay/03.assembly/status/m1000_1/*.tsv"):
    staff_df_m1000 = staff_df_m1000 + [get_status_df(i, "staff_oral")]
for i in glob(staff_oral + "/assay/03.assembly/status/m1000_2/*.tsv"):
    staff_df_m1000 = staff_df_m1000 + [get_status_df(i, "staff_oral")]
staff_asm_df_m1000 = pd.concat(staff_df_m1000)

staff_df_bins_m0 = []
for i in glob(staff_oral + "/assay/05.binning/status/m0/*.tsv"):
    staff_df_bins_m0 = staff_df_bins_m0 + [get_status_df(i, "staff_oral", True)]
staff_asm_bins_df_m0 = pd.concat(staff_df_bins_m0)

#%%
staff_asm_df_m0.head()

#%%
#### pub1 oral
pub1_df_m0 = []
for i in glob(pub1_oral + "/assay/03.assembly/status/m0/*.tsv"):
    pub1_df_m0 = pub1_df_m0 + [get_status_df(i, "pub1_oral")]
pub1_asm_df_m0 = pd.concat(pub1_df_m0)

pub1_df_m1000 = []
for i in glob(pub1_oral + "/assay/03.assembly/status/m1000/*.tsv"):
    pub1_df_m1000 = pub1_df_m1000 + [get_status_df(i, "pub1_oral")]
pub1_asm_df_m1000 = pd.concat(pub1_df_m1000)

pub1_df_bins_m0 = []
for i in glob(pub1_oral + "/assay/05.binning/status/m0/*.tsv"):
    pub1_df_bins_m0 = pub1_df_bins_m0 + [get_status_df(i, "pub1_oral", True)]
pub1_asm_bins_df_m0 = pd.concat(pub1_df_bins_m0)

#%%
#### pub2 oral
pub2_df_m0 = []
for i in glob(pub2_oral + "/assay/03.assembly/status/m0/*.tsv"):
    pub2_df_m0 = pub2_df_m0 + [get_status_df(i, "pub2_oral")]
pub2_asm_df_m0 = pd.concat(pub2_df_m0)

pub2_df_m1000 = []
for i in glob(pub2_oral + "/assay/03.assembly/status/m1000/*.tsv"):
    pub2_df_m1000 = pub2_df_m1000 + [get_status_df(i, "pub2_oral")]
pub2_asm_df_m1000 = pd.concat(pub2_df_m1000)

pub2_df_bins_m0 = []
for i in glob(pub2_oral + "/assay/05.binning/status/m0/*.tsv"):
    pub2_df_bins_m0 = pub2_df_bins_m0 + [get_status_df(i, "pub2_oral", True)]
pub2_asm_bins_df_m0 = pd.concat(pub2_df_bins_m0)

#%%
#### yunnan oral
yunnan_df_m0 = []
for i in glob(yunnan_oral + "/assay/03.assembly/status/m0/*.tsv"):
    yunnan_df_m0 = yunnan_df_m0 + [get_status_df(i, "yunnan_oral")]
yunnan_asm_df_m0 = pd.concat(yunnan_df_m0)

yunnan_df_m1000 = []
for i in glob(yunnan_oral + "/assay/03.assembly/status/m1000/*.tsv"):
    yunnan_df_m1000 = yunnan_df_m1000 + [get_status_df(i, "yunnan_oral")]
yunnan_asm_df_m1000 = pd.concat(yunnan_df_m1000)

yunnan_df_bins_m0 = []
for i in glob(yunnan_oral + "/assay/05.binning/status/m0/*.tsv"):
    yunnan_df_bins_m0 = yunnan_df_bins_m0 + [get_status_df(i, "yunnan_oral", True)]
yunnan_asm_bins_df_m0 = pd.concat(yunnan_df_bins_m0)

#%%
yunnan_asm_bins_df_m0.head()

#%%
[sum(staff_asm_df_m0.duplicated("sample_id")),
 sum(staff_asm_df_m1000.duplicated("sample_id")),
 sum(staff_asm_bins_df_m0.duplicated("sample_id")),
 sum(pub1_asm_df_m0.duplicated("sample_id")),
 sum(pub1_asm_df_m1000.duplicated("sample_id")),
 sum(pub1_asm_bins_df_m0.duplicated("sample_id")),
 sum(pub2_asm_df_m0.duplicated("sample_id")),
 sum(pub2_asm_df_m1000.duplicated("sample_id")),
 sum(pub2_asm_bins_df_m0.duplicated("sample_id")),
 sum(yunnan_asm_df_m0.duplicated("sample_id")),
 sum(yunnan_asm_df_m1000.duplicated("sample_id")),
 sum(yunnan_asm_bins_df_m0.duplicated("sample_id"))]
# why all bins have duplicated?
# the assembly summary program may have problem

#%%
asm_df_m0 = pd.concat([staff_asm_df_m0, pub1_asm_df_m0, pub2_asm_df_m0, yunnan_asm_df_m0],
                      ignore_index=True)\
              .drop_duplicates("sample_id")
asm_df_m1000 = pd.concat([staff_asm_df_m1000, pub1_asm_df_m1000, pub2_asm_df_m1000, yunnan_asm_df_m1000],
                         ignore_index=True)\
                 .drop_duplicates("sample_id")
asm_bins_df_m0 = pd.concat([staff_asm_bins_df_m0, pub1_asm_bins_df_m0, pub2_asm_bins_df_m0, yunnan_asm_bins_df_m0],
                           ignore_index=True)\
                   .drop_duplicates("bin_id")

#%%
[len(asm_df_m0), len(asm_df_m1000), len(asm_bins_df_m0)]
## total 163718 bins

#%%
asm_df_m0_sorted = asm_df_m0.sort_values(by=["scaf_N50", "scaf_bp", "n_scaffolds"], ascending=[False, False, True])
asm_df_m1000_sorted = asm_df_m1000.sort_values(by=["scaf_N50", "scaf_bp", "n_scaffolds"], ascending=[False, False, True])
asm_bins_df_m0_sorted = asm_bins_df_m0.sort_values(by=["scaf_N50", "scaf_bp", "n_scaffolds"], ascending=[False, False, True])

#%%
asm_df_m0_sorted.to_csv(os.path.join(save_dir, "assembly_summary.m0.tsv"), sep='\t', index=False)
asm_df_m1000_sorted.to_csv(os.path.join(save_dir, "assembly_summary.m1000.tsv"), sep='\t', index=False)
asm_bins_df_m0_sorted.to_csv(os.path.join(save_dir, "assembly_summary.bins.m0.tsv"), sep='\t', index=False)

################################################################################################################################
#%%
staff_asm_df_m1000_sorted.to_csv(staff_oral + "/assay/03.assembly/status/asm_status_m1000_sorted_summary.tsv", sep='\t', index=False)

#%%
staff_asm_df_m0.describe().style

#%%
staff_asm_df_m1000.describe().style

#%%
staff_asm_bins_df_m0_dedup.head()

#%%
staff_asm_bins_df_m0_dedup.describe().style

#%%
staff_asm_df_m0.describe().to_csv(staff_oral + "/assay/03.assembly/status/assembly_status_m0.tsv", sep='\t', index=True)

#%%
staff_asm_df_m1000.describe().to_csv(staff_oral + "/assay/03.assembly/status/assembly_status_m1000.tsv", sep='\t', index=True)

#%%
staff_asm_bins_df_m0_dedup.describe().to_csv(staff_oral + "/assay/05.binning/status/assembly_status_bins_m0.tsv", sep='\t', index=True)

#%%
staff_asm_df_m0["sample_site"] = staff_asm_df_m0.apply(sample_site, axis=1)

#%%
staff_asm_df_m1000["sample_site"] = staff_asm_df_m1000.apply(sample_site, axis=1)

#%%
staff_asm_bins_df_m0_dedup["sample_site"] = staff_asm_bins_df_m0_dedup.apply(sample_site, axis=1)

#%%
staff_asm_df_m0_l = staff_asm_df_m0.melt(id_vars=["sample_id", "sample_site"],
                                         var_name="class", value_name = "value")

#%%
staff_asm_df_m1000_l = staff_asm_df_m1000.melt(id_vars=["sample_id", "sample_site"],
                                               var_name="class", value_name = "value")

#%%
staff_asm_bins_df_m0_dedup.head()

#%%
staff_asm_bins_df_m0_l = staff_asm_bins_df_m0_dedup.melt(id_vars=["sample_id", "bin_id", "sample_site"],
                                                         var_name="class", value_name = "value")

#%%
len(staff_asm_bins_df_m0_l)

#%%
plotnine.options.figure_size = (12, 12)
staff_asm_status_plot_m0 = (ggplot(staff_asm_df_m0_l, aes(x='sample_site', y='value'))
                         + geom_boxplot(aes(fill='class', colour='class'))
                         + facet_wrap('~ class', scales='free')
                         + ggtitle("2674 samples assembly statistics metrics summary(m0)(staff_oral)"))
staff_asm_status_plot_m0.save(staff_oral + "/assay/03.assembly/status/assembly_statistics_metrics_summary_2674S_m0.pdf")
staff_asm_status_plot_m0

#%%
plotnine.options.figure_size = (12, 12)
staff_asm_status_plot_m1000 = (ggplot(staff_asm_df_m1000_l,
                                      aes(x='sample_site', y='value'))
                         + geom_boxplot(aes(fill='class', colour='class'))
                         + facet_wrap('~ class', scales='free')
                         + ggtitle("2674 samples assembly statistics metrics summary(m1000)(staff_oral)"))
staff_asm_status_plot_m1000.save(staff_oral + "/assay/03.assembly/status/assembly_statistics_metrics_summary_2674S_m1000.pdf")
staff_asm_status_plot_m1000

#%%
staff_asm_bins_df_m0_l.head()

#%%
plotnine.options.figure_size = (12, 12)
staff_asm_bins_status_plot_m0 = (ggplot(staff_asm_bins_df_m0_l,
                                        aes(x='sample_site', y='value'))
                           + geom_boxplot(aes(fill='class', colour='class'))
                           + facet_wrap('~ class', scales='free')
                           + ggtitle("126872 bins assembly statistics metrics summary(m1500)(staff_oral)"))
staff_asm_bins_status_plot_m0.save(staff_oral + "/assay/05.binning/status/bins_assembly_statistics_metrics_summary_126872_m1500.pdf")
staff_asm_bins_status_plot_m0
#%%
staff_checkm_summary = pd\
    .read_csv(staff_oral + "/assay/06.checkm/status/checkm_ccsh_summary.tsv",
              sep='\t',
              nrows=2675)

#%%
staff_checkm_out_merge = pd\
    .read_csv(staff_oral + "/assay/06.checkm/status/checkm_ccsh_out.txt",
              sep='\t')
staff_checkm_summary.describe()\
    .to_csv(staff_oral + "/assay/06.checkm/status/checkm_ccsh_summary_describe.tsv",
            sep='\t')
staff_checkm_out_merge.describe()\
    .to_csv(staff_oral + "/assay/06.checkm/status/checkm_ccsh_out_describe.tsv",
            sep='\t')
staff_checkm_summary["sample_site"] = staff_checkm_summary.apply(sample_site, axis=1)
staff_checkm_out_merge["sample_site"] = staff_checkm_out_merge.apply(sample_site2, axis=1)
staff_checkm_summary_l = staff_checkm_summary.melt(id_vars=["sample_id", "sample_site"],
                                                   var_name="class", value_name="value")
staff_checkm_out_merge_l = (staff_checkm_out_merge.loc[:, ["bin_id", "sample_site",
                                                           "completeness",
                                                           "contamination",
                                                           "strain_heterogeneity"]]
                                                  .melt(id_vars=["bin_id", "sample_site"],
                                                        var_name="class", value_name="value"))
plotnine.options.figure_size = (20, 20)
staff_checkm_summary_plot = (ggplot(staff_checkm_summary_l, aes(x='sample_site', y='value'))
                          + geom_boxplot(aes(fill='class', colour='class'))
                          + facet_wrap('~ class', scales='free')
                          + ggtitle("2675 samples checkm statistics metrics summary(staff_oral)"))
staff_checkm_summary_plot.save(staff_oral + "/assay/06.checkm/status/checkm_statistics_metrics_summary_2675S.pdf")
staff_checkm_summary_plot
plotnine.options.figure_size = (8, 6)
staff_checkm_out_merge_plot = (ggplot(staff_checkm_out_merge_l, aes(x='sample_site', y='value'))
                            + geom_boxplot(aes(fill='class', colour='class'))
                            + facet_wrap('~ class', scales='free')
                            + ggtitle("126872 bins statistics metrics summary(staff oral)"))
staff_checkm_out_merge_plot.save(staff_oral + "/assay/06.checkm/status/checkm_statistics_metrics_summary_126872.pdf")
staff_checkm_out_merge_plot
staff_checkm_out_merge_high_quality = staff_checkm_out_merge.query('completeness > 90 & contamination < 5 & strain_heterogeneity < 0.5')
staff_checkm_out_merge_medium_quality = staff_checkm_out_merge.query('completeness <= 90 & completeness > 50 & contamination < 5')
staff_checkm_out_merge_high_medium_df = pd.concat([staff_checkm_out_merge_high_quality, staff_checkm_out_merge_medium_quality])
staff_asm_bins_df_m0_high_quality = staff_asm_bins_df_m0_dedup[staff_asm_bins_df_m0_dedup.bin_id.isin(staff_checkm_out_merge_high_quality.bin_id)]
staff_asm_bins_df_m0_high_medium_quality = staff_asm_bins_df_m0_dedup[staff_asm_bins_df_m0_dedup.bin_id.isin(staff_checkm_out_merge_high_medium_df.bin_id)]
yunnan_asm_df_m0_sorted = yunnan_asm_df_m0.sort_values(by=["scaf_N50", "scaf_bp", "n_scaffolds"], ascending=[False, False, True])
yunnan_asm_df_m0_sorted.to_csv(yunnan_oral + "/assay/03.assembly/status/asm_status_m0_sorted_summary.tsv", sep='\t', index=False)
yunnan_asm_df_m1000_sorted = yunnan_asm_df_m1000.sort_values(by=["scaf_N50", "scaf_bp", "n_scaffolds"], ascending=[False, False, True])
yunnan_asm_df_m1000_sorted.to_csv(yunnan_oral + "/assay/03.assembly/status/asm_status_m1000_sorted_summary.tsv", sep='\t', index=False)
yunnan_asm_df_m0.describe().to_csv(yunnan_oral + "/assay/03.assembly/status/assembly_status_m0.tsv", sep='\t', index=True)
yunnan_asm_df_m1000.describe().to_csv(yunnan_oral + "/assay/03.assembly/status/assembly_status_m1000.tsv", sep='\t', index=True)
yunnan_asm_bins_df_m0.describe().to_csv(yunnan_oral + "/assay/05.binning/status/assembly_status_bins_m0.tsv", sep='\t', index=True)
yunnan_asm_df_m0["sample_site"] = yunnan_asm_df_m0.apply(sample_site, axis=1)
yunnan_asm_df_m1000["sample_site"] = yunnan_asm_df_m1000.apply(sample_site, axis=1)
yunnan_asm_bins_df_m0["sample_site"] =yunnan_asm_bins_df_m0.apply(sample_site, axis=1)
yunnan_asm_df_m0_l = yunnan_asm_df_m0.melt(id_vars=["sample_id", "sample_site"],
                                           var_name="class", value_name = "value")
yunnan_asm_df_m1000_l = yunnan_asm_df_m1000.melt(id_vars=["sample_id", "sample_site"],
                                                 var_name="class", value_name = "value")
yunnan_asm_bins_df_m0_l = yunnan_asm_bins_df_m0.melt(id_vars=["sample_id", "bin_id", "sample_site"],
                                                     var_name="class", value_name = "value")
plotnine.options.figure_size = (12, 12)
yunnan_asm_status_plot_m0 = (ggplot(yunnan_asm_df_m0_l, aes(x='sample_site', y='value'))
                         + geom_boxplot(aes(fill='class', colour='class'))
                         + facet_wrap('~ class', scales='free')
                         + ggtitle("671 samples assembly statistics metrics summary(m0)(yunnan_oral)"))
yunnan_asm_status_plot_m0.save(yunnan_oral + "/assay/03.assembly/status/assembly_statistics_metrics_summary_671S_m0.pdf")
yunnan_asm_status_plot_m0
plotnine.options.figure_size = (12, 12)
yunnan_asm_status_plot_m1000 = (ggplot(yunnan_asm_df_m1000_l, aes(x='sample_site', y='value'))
                         + geom_boxplot(aes(fill='class', colour='class'))
                         + facet_wrap('~ class', scales='free')
                         + ggtitle("671 samples assembly statistics metrics summary(m1000)(yunnan_oral)"))
yunnan_asm_status_plot_m1000.save(yunnan_oral + "/assay/03.assembly/status/assembly_statistics_metrics_summary_671S_m1000.pdf")
yunnan_asm_status_plot_m1000
plotnine.options.figure_size = (12, 12)
yunnan_asm_bins_status_plot_m0 = (ggplot(yunnan_asm_bins_df_m0_l, aes(x='sample_site', y='value'))
                           + geom_boxplot(aes(fill='class', colour='class'))
                           + facet_wrap('~ class', scales='free')
                           + ggtitle("12922 bins assembly statistics metrics summary(m1500)(yunnan_oral)"))
yunnan_asm_bins_status_plot_m0.save(yunnan_oral + "/assay/05.binning/status/bins_assembly_statistics_metrics_summary_12922_m1500.pdf")
yunnan_asm_bins_status_plot_m0
yunnan_oral + "/assay/06.checkm/status/checkm_ccsh_summary.tsv"
get_ipython().system('wc -l /hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/yunnan_oral/assay/06.checkm/status/checkm_ccsh_summary.tsv')
yunnan_checkm_summary = pd.read_csv(yunnan_oral + "/assay/06.checkm/status/checkm_ccsh_summary.tsv", sep='\t', nrows=671)
yunnan_checkm_out_merge = pd.read_csv(yunnan_oral + "/assay/06.checkm/status/checkm_ccsh_out.txt", sep='\t')
yunnan_checkm_summary.describe().to_csv(yunnan_oral + "/assay/06.checkm/status/checkm_ccsh_summary_describe.tsv", sep='\t')
yunnan_checkm_out_merge.describe().to_csv(yunnan_oral + "/assay/06.checkm/status/checkm_ccsh_out_describe.tsv", sep='\t')
yunnan_checkm_summary["sample_site"] = yunnan_checkm_summary.apply(sample_site, axis=1)
yunnan_checkm_out_merge["sample_site"] =yunnan_checkm_out_merge.apply(sample_site2, axis=1)
yunnan_checkm_summary_l = yunnan_checkm_summary.melt(id_vars=["sample_id", "sample_site"],
                                                   var_name="class", value_name="value")
yunnan_checkm_out_merge_l = (yunnan_checkm_out_merge.loc[:, ["bin_id", "sample_site",
                                                           "completeness",
                                                           "contamination",
                                                           "strain_heterogeneity"]]
                                                  .melt(id_vars=["bin_id", "sample_site"],
                                                        var_name="class", value_name="value"))
plotnine.options.figure_size = (20, 20)
yunnan_checkm_summary_plot = (ggplot(yunnan_checkm_summary_l, aes(x='sample_site', y='value'))
                          + geom_boxplot(aes(fill='class', colour='class'))
                          + facet_wrap('~ class', scales='free')
                          + ggtitle("671 samples checkm statistics metrics summary(yunnan_oral)"))
yunnan_checkm_summary_plot.save(yunnan_oral + "/assay/06.checkm/status/checkm_statistics_metrics_summary_671S.pdf")
yunnan_checkm_summary_plot
plotnine.options.figure_size = (8, 6)
yunnan_checkm_out_merge_plot = (ggplot(yunnan_checkm_out_merge_l, aes(x='sample_site', y='value'))
                            + geom_boxplot(aes(fill='class', colour='class'))
                            + facet_wrap('~ class', scales='free')
                            + ggtitle("12922 bins statistics metrics summary(yunnan oral)"))
yunnan_checkm_out_merge_plot.save(staff_oral + "/assay/06.checkm/status/checkm_statistics_metrics_summary_12922.pdf")
yunnan_checkm_out_merge_plot
yunnan_checkm_out_merge_high_quality = yunnan_checkm_out_merge.query('completeness > 90 & contamination < 5 & strain_heterogeneity < 0.5')
yunnan_checkm_out_merge_medium_quality = yunnan_checkm_out_merge.query('completeness <= 90 & completeness > 50 & contamination < 5')
yunnan_checkm_out_merge_high_medium_df = pd.concat([yunnan_checkm_out_merge_high_quality, yunnan_checkm_out_merge_medium_quality])
yunnan_asm_bins_df_m0_high_quality = yunnan_asm_bins_df_m0[yunnan_asm_bins_df_m0.bin_id.isin(yunnan_checkm_out_merge_high_quality.bin_id)]
yunnan_asm_bins_df_m0_high_medium_quality = yunnan_asm_bins_df_m0[yunnan_asm_bins_df_m0.bin_id.isin(yunnan_checkm_out_merge_high_medium_df.bin_id)]
pub1_oral + "/assay/06.checkm/status/checkm_ccsh_summary.tsv"
get_ipython().system('wc -l /hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/pub_oral/assay/06.checkm/status/checkm_ccsh_summary.tsv')


#%%
pub1_checkm_summary = pd.read_csv(pub1_oral + "/assay/06.checkm/status/checkm_ccsh_summary.tsv", sep='\t', nrows=508)


#%%
pub1_checkm_out_merge = pd.read_csv(pub1_oral + "/assay/06.checkm/status/checkm_ccsh_out.txt", sep='\t')


#%%
pub1_checkm_out_merge_high_quality = pub1_checkm_out_merge.query('completeness > 90 &                                                                       contamination < 5 &                                                                       strain_heterogeneity < 0.5')


#%%
len(pub1_checkm_out_merge_high_quality)


#%%
pub1_checkm_out_merge_medium_quality = pub1_checkm_out_merge.query('completeness <= 90 &                                                                        completeness > 50 &                                                                        contamination < 5')


#%%
len(pub1_checkm_out_merge_medium_quality)


#%%
pub1_checkm_out_merge_high_medium_df = pd.concat([pub1_checkm_out_merge_high_quality, pub1_checkm_out_merge_medium_quality])


#%%
pub1_asm_bins_df_m0_high_quality = pub1_asm_bins_df_m0[pub1_asm_bins_df_m0.bin_id.isin(pub1_checkm_out_merge_high_quality.bin_id)]


#%%
pub1_asm_bins_df_m0_high_medium_quality = pub1_asm_bins_df_m0[pub1_asm_bins_df_m0.bin_id.isin(pub1_checkm_out_merge_high_medium_df.bin_id)]


#%%
pub1_asm_bins_df_m0_high_medium_quality.describe().style


#%%
###### pub 2 oral (ra_oral)


#%%



#%%
len(pub2_asm_bins_df_m0)


#%%
pub2_oral + "/assay/06.checkm/status/checkm_ccsh_summary.tsv"


#%%
get_ipython().system('wc -l /hwfssz1/ST_META/F16ZQSB1SY2794_RO_GRJ/RO2/assay/06.checkm/status/checkm_ccsh_summary.tsv')


#%%
pub2_checkm_summary = pd.read_csv(pub2_oral + "/assay/06.checkm/status/checkm_ccsh_summary.tsv", sep='\t', nrows=294)


#%%
pub2_checkm_out_merge = pd.read_csv(pub2_oral + "/assay/06.checkm/status/checkm_ccsh_out.txt", sep='\t')


#%%
pub2_checkm_out_merge_high_quality = pub2_checkm_out_merge.query('completeness > 90 &                                                                       contamination < 5 &                                                                       strain_heterogeneity < 0.5')


#%%
len(pub2_checkm_out_merge_high_quality)


#%%
pub2_checkm_out_merge_medium_quality = pub2_checkm_out_merge.query('completeness <= 90 &                                                                        completeness > 50 &                                                                        contamination < 5')


#%%
len(pub2_checkm_out_merge_medium_quality)


#%%
pub2_checkm_out_merge_high_medium_df = pd.concat([pub2_checkm_out_merge_high_quality, pub2_checkm_out_merge_medium_quality])


#%%
pub2_asm_bins_df_m0_high_quality = pub2_asm_bins_df_m0[pub2_asm_bins_df_m0.bin_id.isin(pub2_checkm_out_merge_high_quality.bin_id)]


#%%
len(pub2_asm_bins_df_m0_high_quality)


#%%
pub2_asm_bins_df_m0_high_medium_quality = pub2_asm_bins_df_m0[pub2_asm_bins_df_m0.bin_id.isin(pub2_checkm_out_merge_high_medium_df.bin_id)]


#%%
pub2_asm_bins_df_m0_high_medium_quality.head()


#%%
pub2_asm_bins_df_m0_high_medium_quality.describe().style


#%%
############################################################################################################################################


#%%
# merge analysis


#%%
sum(staff_asm_bins_df_m0.duplicated("bin_id"))

staff_asm_bins_df_m0_dedup
yunnan_asm_bins_df_m0
pub1_asm_bins_df_m0
pub2_asm_bins_df_m0
#%%
staff_asm_bins_df_m0_dedup["project"] = "staff_oral"
yunnan_asm_bins_df_m0["project"] = "yunnan_oral" 
pub1_asm_bins_df_m0["project"] = "pub1_oral"
pub2_asm_bins_df_m0["project"] = "pub2_oral"


#%%
oral_asm_bins_df_m0 = pd.concat([staff_asm_bins_df_m0_dedup, yunnan_asm_bins_df_m0, pub1_asm_bins_df_m0, pub2_asm_bins_df_m0])


#%%
len(oral_asm_bins_df_m0)


#%%
oral_asm_bins_df_m0.head()


#%%
oral_asm_bins_df_m0_l = oral_asm_bins_df_m0.loc[:, ["bin_id", "sample_id", "n_scaffolds", "project", "scaf_N50", "scaf_bp"]]                                           .melt(id_vars=["sample_id", "bin_id", "project"],
                                                 var_name="class", value_name = "value")


#%%
len(oral_asm_bins_df_m0_l)


#%%
oral_asm_bins_df_m0_l.head()


#%%
508 + 2675 + 294 + 671


#%%
# all bins
plotnine.options.figure_size = (10, 6)
oral_asm_status_plot_all_bins = (ggplot(oral_asm_bins_df_m0_l, aes(x='project', y='value'))
                                + geom_boxplot(aes(fill='class', colour='class'))
                                + facet_wrap('~ class', scales='free')
                                + ggtitle("4148 oral samples assembly statistics metrics summary of 163718 bins")
                                + theme(axis_text_x=element_text(angle=90)))
oral_asm_status_plot_all_bins.save("4148_oral_samples_assembly_statistics_metrics_summary_all_bins.pdf")
oral_asm_status_plot_all_bins


#%%
# high quality and medium bins
oral_asm_bins_df_m0_high_quality = pd.concat([staff_asm_bins_df_m0_high_quality,
                                              yunnan_asm_bins_df_m0_high_quality,
                                              pub1_asm_bins_df_m0_high_quality,
                                              pub2_asm_bins_df_m0_high_quality])
oral_asm_bins_df_m0_high_medium_quality = pd.concat([staff_asm_bins_df_m0_high_medium_quality,
                                              yunnan_asm_bins_df_m0_high_medium_quality,
                                              pub1_asm_bins_df_m0_high_medium_quality,
                                              pub2_asm_bins_df_m0_high_medium_quality])


#%%
len(oral_asm_bins_df_m0_high_quality)


#%%
len(oral_asm_bins_df_m0_high_medium_quality)


#%%
oral_asm_bins_df_m0_high_quality_l = oral_asm_bins_df_m0_high_quality.loc[:, ["bin_id", "sample_id", "n_scaffolds", "project", "scaf_N50", "scaf_bp"]]                                           .melt(id_vars=["sample_id", "bin_id", "project"],
                                                 var_name="class", value_name = "value")
oral_asm_bins_df_m0_high_medium_quality_l = oral_asm_bins_df_m0_high_medium_quality.loc[:, ["bin_id", "sample_id", "n_scaffolds", "project", "scaf_N50", "scaf_bp"]]                                           .melt(id_vars=["sample_id", "bin_id", "project"],
                                                 var_name="class", value_name = "value")


#%%
# high quality bins
plotnine.options.figure_size = (10, 6)
oral_asm_status_plot_high_quality_bins = (ggplot(oral_asm_bins_df_m0_high_quality_l, aes(x='project', y='value'))
                                + geom_boxplot(aes(fill='class', colour='class'))
                                + facet_wrap('~ class', scales='free')
                                + ggtitle("4148 oral samples assembly statistics metrics summary of 5145 high quality bins")
                                + theme(axis_text_x=element_text(angle=90)))
oral_asm_status_plot_high_quality_bins.save("4148_oral_samples_assembly_statistics_metrics_summary_5145_high_quality.pdf")
oral_asm_status_plot_high_quality_bins


#%%
oral_asm_bins_df_m0_high_medium_quality_l.head()


#%%
# high quality and medium bins
plotnine.options.figure_size = (10, 6)
oral_asm_status_plot_high_medium_quality_bins = (ggplot(oral_asm_bins_df_m0_high_medium_quality_l, aes(x='project', y='value'))
                                                + geom_boxplot(aes(fill='class', colour='class'))
                                                + facet_wrap('~ class', scales='free')
                                                + ggtitle("4148 oral samples assembly statistics metrics summary of 35895 high and medium quality bins")
                                                + theme(axis_text_x=element_text(angle=90)))
oral_asm_status_plot_high_medium_quality_bins.save("4148_oral_samples_assembly_statistics_metrics_summary_35895_high_and_medium_quality.pdf")
oral_asm_status_plot_high_medium_quality_bins


#%%
oral_asm_bins_df_m0_high_quality.describe()


#%%
oral_asm_bins_df_m0_high_medium_quality.describe()


#%%
oral_asm_bins_df_m0_high_quality.project.value_counts()


#%%
3973 / 5145


#%%
oral_asm_bins_df_m0_high_medium_quality.project.value_counts()


#%%
26789 / 35895


#%%
# mapping analysis


#%%
[staff_oral, yunnan_oral, pub1_oral, pub2_oral]


#%%
get_ipython().system('ls -l /hwfssz5/ST_META/P18Z10200N0127_MA/zhujie/sgb_oral/assay/04.alignment/status/flagstat.txt')


#%%
staff_flagstat = staff_oral + "/assay/04.alignment/status/flagstat.txt"
yunnan_flagstat = yunnan_oral + "/assay/04.alignment/status/flagstat.txt"
pub1_flagstat = pub1_oral + "/assay/04.alignment/status/flagstat.txt"
pub2_flagstat = pub2_oral + "/assay/04.alignment/status/flagstat.txt"


#%%
pub2_flagstat


#%%
get_ipython().system('head /hwfssz1/ST_META/F16ZQSB1SY2794_RO_GRJ/RO2/assay/04.alignment/status/flagstat.txt')


#%%
def alignment_status(flagstat, project_):
    df = (pd.read_csv(flagstat, sep='\t')
            .assign(project=project_)
            .loc[:, ['sample_id', 'project', 'mapped_rate', 'paired_rate', 'singletons_rate']]
            .sort_values('mapped_rate')
            .melt(id_vars=['sample_id', 'project'], var_name='mapping', value_name='rate'))
    return df


#%%
staff_flagstat_df = alignment_status(staff_flagstat, "staff_oral")
yunnan_flagstat_df = alignment_status(yunnan_flagstat, "yunnan_oral")
pub1_flagstat_df = alignment_status(pub1_flagstat, "pub1_oral")
pub2_flagstat_df = alignment_status(pub2_flagstat, "pub2_oral")


#%%
oral_flagstat_df = pd.concat([staff_flagstat_df, yunnan_flagstat_df, pub1_flagstat_df, pub2_flagstat_df])


#%%
oral_flagstat_df.head()


#%%
len(oral_flagstat_df)


#%%
len(oral_flagstat_df.query('mapping == "mapped_rate" and project == "pub2_oral"'))


#%%
oral_flagstat_p = (ggplot(oral_flagstat_df.query('mapping == "mapped_rate" and project == "pub2_oral"'),
                          aes(x='sample_id', y='rate'))
                   + geom_point(aes(color='mapping', fill='mapping'))
                   + facet_wrap('~ project')
                   + theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank()))
                   
oral_flagstat_p


#%%
# high quality don't think strain heretogenity

# high quality think strain heretogenity
