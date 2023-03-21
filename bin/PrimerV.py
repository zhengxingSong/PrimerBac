import pandas as pd
import os
import subprocess
import argparse
import re

def configP(dat):
    with open(dat,"r") as f:
        input_str = f.read()
    input_dict = {}
    for line in input_str.strip().split('\n'):
        key, value = line.split(':')
        try:
            input_dict[key] = int(value.strip())
        except ValueError:
            input_dict[key] = value.strip()
    return(input_dict)

def flagTrans(Num):
    flagL = []
    NumL = list(reversed("{0:b}".format(Num)))
    for i in range(len(NumL)):
        if NumL[i] != "0":
            flagL.append(2**i)
    return(flagL)

def bowtieR(Infile,id2sp):
    with open(Infile,"r") as f:
        data = f.readlines()

    data = [x.split("\t") for x in data if not x.startswith("@")]
    dat = pd.DataFrame(data)

    flag_need = [x for x in dat[1].unique().tolist() if 4 not in flagTrans(int(x))]

    dat = dat[dat[1].isin(flag_need)]
    dat = dat.iloc[:, [0,2,8]]
    dat.columns = ["Primer", "ref", "primerL"]
    dat["primerL"] = dat["primerL"].str.extract("(\d+)").astype(int)
    dat["Orient"] = dat["Primer"].str.split("_").str[-1]
    dat["primer"] = dat["Primer"].str.replace("_[LR]$", "", regex=True)
    dat["ref"] = dat["ref"].apply(lambda x:x.split("_")[0] if x.startswith("MGYG") else x)
    dat["Ref_sp"] = dat["ref"].map(id2sp).fillna("UnKnown")

    dat_next = dat.groupby("primer").agg(
        primerL=("primerL", "first"),
        Orient=("Orient", lambda x: ",".join(sorted(set(x)))),
        Ref_sp=("Ref_sp", lambda x: ",".join(sorted(set(x)))),
        ref=("ref", lambda x: ",".join(sorted(set(x)))),
    ).reset_index()

    dat_next = dat_next[dat_next["primerL"] != 0]
    return(dat_next)

def SpeciesR(Infile,id2sp,bac,BacNum):
    dat = bowtieR(Infile,id2sp)
    dat["Target"] = dat["Ref_sp"].apply(lambda x:len([y for y in x.split(",") if re.findall(bac,y)])/len(x.split(","))*100)
    dat["Senstive"] = dat["Ref_sp"].apply(lambda x:len([y for y in x.split(",") if re.findall(bac,y)])/BacNum[bac]*100)
    dat = dat.drop(["Orient","primerL","ref"],axis=1)
    dat = dat.query("Target!=0 & Senstive!=0")
    return(dat)

class Bowtie:
    def __init__(self, index_path, reads_path, output_path,bac):
        self.index_path = index_path
        self.reads_path = reads_path
        self.output_path = output_path
        self.bac = bac

    def run(self):
        cmd = ['bowtie2', '-p', '4', '-f',"-a","--seed","123","-x",self.index_path,
                '-1',f"{self.reads_path}/{self.bac}_1.fa","-2",f"{self.reads_path}/{self.bac}_2.fa",
                "-S",f"{self.output_path}/{self.bac}.sam"]
        subprocess.call(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        required=True,
        help="The input yaml file.",
    )
    args = parser.parse_args()
    conF = args.config
    config = configP(conF)
    default = configP(os.path.join(os.path.dirname(os.path.abspath(__file__)),"default.yaml"))
    df_intense = pd.read_table(default["Taxanno"])
    id2sp = dict(zip(df_intense["Gene_ID"],df_intense["Species"]))
    Ge2Spnum =df_intense[["Genus","Species"]].drop_duplicates()["Genus"].value_counts().to_dict()
    Sp2StNum = df_intense["Species"].value_counts().to_dict()
    BacNum = Ge2Spnum if config["Taxonomy"] == "Genus" else Sp2StNum
    bowtieIn = default["BowIndex"]
    bac = config["Bacteria"].replace(" ","_")
    df_Meta = pd.read_table(os.path.join(config["Result"],bac,"%s_metadata.txt"%bac))
    for index,value in df_Meta.iterrows():
        InPath = os.path.join(config["Result"],value["Genus"],value["Species"].replace(" ","_"))
        # bt = Bowtie(bowtieIn,os.path.join(InPath,"Working_Direct"),InPath,bac)
        # bt.run()
        dat = SpeciesR(os.path.join(InPath,f"{bac}.sam"),id2sp,config["Bacteria"],BacNum)
        dat.to_csv(os.path.join(InPath,"FinalP.txt"),index=False,sep="\t")

if __name__ == "__main__":
    main()