import pandas as pd

NON_MOD_BASES = ["pAr", "pUr", "pGr", "pCr"]

def find_mod_pos(row):
    oligo_seq_list = ("p"+row).split("-")[:-1]
    return [(x,i) for i, x in enumerate(oligo_seq_list) if x not in NON_MOD_BASES]

def generate_bedrmod(filename, filename_mapping):
    df = pd.read_csv(filename)
    df_mapping = pd.read_csv(filename_mapping, sep = "\t")
    
    df = df[(df["Best ASR"] <= 1) & (df["Conf. Score"] >= 100)]
    
    df["pos_start"] = df["Positions"].apply(lambda x: int(x.split("-")[0]))
    df["pos_end"] = df["Positions"].apply(lambda x: int(x.split("-")[1]))
    
    df["mod"] = df["Oligo Sequence"].apply(find_mod_pos)
    
    df2 = df[df["mod"].map(len)>0][["pos_start", "mod", "Conf. Score"]]
    
    df2 = df2.explode("mod")
    df2["mod_name"] = df2["mod"].apply(lambda x: x[0])
    df2["mod_idx"] = df2["mod"].apply(lambda x: x[1])
    
    df2["chromEnd"] = df2["pos_start"]+df2["mod_idx"]
    df2["chromStart"] = df2["chromEnd"]-1
    
    df2 = df2.drop(["mod", "pos_start", "mod_idx"], axis = 1)
    
    df3 = df2.groupby(["mod_name", 
                       "chromStart", 
                       "chromEnd"]).agg(coverage = ("mod_name", "size"),
                                        score = ("Conf. Score", "mean")).reset_index()
    df3["frequency"] = (df3["coverage"]/df.shape[0])*100
    
    df3 = pd.merge(df3, df_mapping, how = "left", left_on = "chromEnd", right_on = "pos")
    df3 = df3.drop(["mod_name", "pos"], axis = 1)
    
    df3["thickStart"] = df3["chromStart"]
    df3["thickEnd"] = df3["chromEnd"]
    
    df3["chrom"] = filename.split("_")[2]
    df3["strand"] = "."
    df3["itemRgb"] = 0
       
    df_final = df3.loc[:, ["chrom",
                           "chromStart",
                           "chromEnd",
                           "name",
                           "score",
                           "strand",
                           "thickStart",
                           "thickEnd",
                           "itemRgb",
                           "coverage",
                           "frequency"]
                       ].sort_values(by = ["chromStart"]).reset_index(drop = True)
    
    return df_final

files = [
         "bpf_mapping/HRPC_5.8S_rRNA.csv",
         "bpf_mapping/HRPC_18S_rRNA_02.csv",
         "bpf_mapping/HRPC_28S_rRNA.csv"
         ]

mapping_files = [
                "bpf_mapping/update_meta_5.8S.tsv",
                "bpf_mapping/update_meta_18S.tsv",
                "bpf_mapping/update_meta_28S.tsv"
                ]

dfs = [generate_bedrmod(files[i], mapping_files[i]) for i in range(len(files))]

df_final = pd.concat(dfs).reset_index(drop = True)

df_final.to_csv("hrpc_bedrmod.tsv", index = False, header = False, sep = "\t")