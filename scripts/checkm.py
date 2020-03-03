#%%
import os
try:
	os.chdir(os.path.join(os.getcwd(), '/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/igo/assay/checkm'))
	print(os.getcwd())
except:
	pass

#%%
print(os.getcwd())

#%%
import pandas as pd

#%%
pub1_oral_checkm = pd\
	.read_csv("/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/pub_oral/assay/06.checkm/status/checkm_ccsh_out.txt",
			  sep='\t')\
	.assign(project="pub1_oral")
pub2_oral_checkm = pd\
	.read_csv("/hwfssz1/ST_META/F16ZQSB1SY2794_RO_GRJ/RO2/assay/06.checkm/status/checkm_ccsh_out.txt",
			  sep='\t')\
	.assign(project="pub2_oral")
staff_oral_checkm = pd\
	.read_csv("/hwfssz5/ST_META/P18Z10200N0127_MA/zhujie/sgb_oral/assay/06.checkm/status/checkm_ccsh_out.txt",
			  sep='\t')\
	.assign(project="staff_oral")
yunnan_oral_checkm = pd\
	.read_csv("/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/yunnan_oral/assay/06.checkm/status/checkm_ccsh_out.txt",
		      sep='\t')\
	.assign(project="yunnan_oral")

#%%
def parse_bins_path(bins_list):
	df = pd.read_csv(bins_list, sep='\t', names=["bins_fna_path"])
	df = df.assign(bin_id = df.bins_fna_path.apply(lambda x: os.path.basename(x).rpartition(".")[0]))
	return df

#%%
pub1_bins = parse_bins_path("/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/pub_oral/assay/05.binning/status/bins_fa.list")
pub2_bins = parse_bins_path("/hwfssz1/ST_META/F16ZQSB1SY2794_RO_GRJ/RO2/assay/05.binning/status/bins_fa.list")
staff_bins = parse_bins_path("/hwfssz5/ST_META/P18Z10200N0127_MA/zhujie/sgb_oral/assay/05.binning/bins_126872.list")
yunnan_bins = parse_bins_path("/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/yunnan_oral/assay/05.binning/status/bins_fa.list")

#%%
pub1_bins.head()

#%%
bins_out = pd.concat([pub1_bins, pub2_bins, staff_bins, yunnan_bins], ignore_index=True)
checkm_out = pd.concat([pub1_oral_checkm,
                        pub2_oral_checkm,
						staff_oral_checkm,
						yunnan_oral_checkm], ignore_index=True)
checkm_out = pd.merge(checkm_out, bins_out)

#%%
checkm_out.head()

#%%
def mimag_quality_level(row):
	# blueprint using this
	if (row["completeness"] > 90) and (row["contamination"] < 5):
		return "high_quality"
	elif (row["completeness"] > 50) and (row["contamination"] < 10):
		return "medium_quality"
	else:
		return "low_quality"

def sgb_quality_level(row):
	if (row["strain_heterogeneity"] < 0.5) and \
	   (row["completeness"] > 90) and (row["contamination"] < 5):
		return "high_quality"
	elif (row["completeness"] > 50) and (row["contamination"] < 5):
		return "medium_quality"
	else:
		return "low_quality"

def quality_score(row):
	# blueprint using this
	return row["completeness"] - 5 * row["contamination"]


#%%
checkm_out["mimag_quality_level"] = checkm_out.apply(mimag_quality_level, axis=1)
checkm_out["sgb_quality_level"] = checkm_out.apply(sgb_quality_level, axis=1)
checkm_out["quality_score"] = checkm_out.apply(quality_score, axis=1)

#%%
checkm_out.head()

#%%
asm_bins_df = pd.read_csv("/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/igo/assay/assembly/results/assembly_summary.bins.m0.tsv",
                          sep='\t')

#%%
checkm_out = checkm_out.merge(asm_bins_df)

#%%
len(checkm_out)

#%%
checkm_out.to_csv("/hwfssz1/ST_META/P18Z10200N0127_MA/zhujie/igo/assay/checkm/checkm_out.tsv",
                  sep='\t', index=False)

#%%
[len(checkm_out.query('mimag_quality_level == "high_quality"')),
 len(checkm_out.query('mimag_quality_level == "medium_quality"'))]

#%%
[len(checkm_out.query('mimag_quality_level == "high_quality" and quality_score > 50')),
 len(checkm_out.query('mimag_quality_level == "medium_quality" and quality_score > 50'))]

#%%
[len(checkm_out.query('sgb_quality_level == "high_quality"')),
 len(checkm_out.query('sgb_quality_level == "medium_quality"'))]

#%%
[len(checkm_out.query('sgb_quality_level == "high_quality" and quality_score > 50')),
 len(checkm_out.query('sgb_quality_level == "medium_quality" and quality_score > 50'))]

#%%
bins_mimag_high_quality = checkm_out\
	.query('mimag_quality_level == "high_quality"')\
	.loc[:, ["bins_fna_path"]]
bins_mimag_high_quality.to_csv("bins_mimag_high_quality.fna.list",
						 sep='\t', index=False, header=False)

#%%
bins_mimag_medium_quality = checkm_out\
	.query('mimag_quality_level == "medium_quality"')\
	.loc[:, ["bins_fna_path"]]
bins_mimag_medium_quality.to_csv("bins_mimag_medium_quality.fna.list",
                           sep='\t', index=False, header=False)

#%%
bins_mimag_low_quality = checkm_out\
	.query('mimag_quality_level == "low_quality"')\
	.loc[:, ["bins_fna_path"]]
bins_mimag_low_quality.to_csv("bins_mimag_low_quality.fna.list",
                        sep='\t', index=False, header=False)

#%%
bins_mimag_high_medium_quality = pd.concat([bins_mimag_high_quality, bins_mimag_medium_quality],
                                     ignore_index=False)
bins_mimag_medium_low_quality = pd.concat([bins_mimag_medium_quality, bins_mimag_low_quality],
                                    ignore_index=False)

#%%
bins_mimag_high_medium_quality.to_csv("bins_mimag_high_medium_quality.fna.list",
                                sep='\t', index=False, header=False)
bins_mimag_medium_low_quality.to_csv("bins_mimag_medium_low_quality.fna.list",
                                sep='\t', index=False, header=False)

#%%
bins_sgb_high_quality = checkm_out\
	.query('sgb_quality_level == "high_quality"')\
	.loc[:, ["bins_fna_path"]]
bins_sgb_high_quality.to_csv("bins_sgb_high_quality.fna.list",
						 sep='\t', index=False, header=False)

#%%
bins_sgb_medium_quality = checkm_out\
	.query('sgb_quality_level == "medium_quality"')\
	.loc[:, ["bins_fna_path"]]
bins_sgb_medium_quality.to_csv("bins_sgb_medium_quality.fna.list",
                           sep='\t', index=False, header=False)

#%%
bins_sgb_low_quality = checkm_out\
	.query('sgb_quality_level == "low_quality"')\
	.loc[:, ["bins_fna_path"]]
bins_sgb_low_quality.to_csv("bins_sgb_low_quality.fna.list",
                        sep='\t', index=False, header=False)

#%%
bins_sgb_high_medium_quality = pd.concat([bins_sgb_high_quality, bins_sgb_medium_quality],
                                     ignore_index=False)
bins_sgb_medium_low_quality = pd.concat([bins_sgb_medium_quality, bins_sgb_low_quality],
                                    ignore_index=False)

#%%
bins_sgb_high_medium_quality.to_csv("bins_sgb_high_medium_quality.fna.list",
                                sep='\t', index=False, header=False)
bins_sgb_medium_low_quality.to_csv("bins_sgb_medium_low_quality.fna.list",
                                sep='\t', index=False, header=False)
