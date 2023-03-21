import pandas as pd
import numpy as np
import os
import re
import primer3
from joblib import Parallel,delayed
import argparse
import logging

# 设置日志记录器
logging.basicConfig(filename='PrimerD.log', level=logging.INFO, 
                    format='%(asctime)s %(levelname)s: %(message)s')


def Fa2D(InFile):
    FaD = dict()
    with open(InFile,"r") as f:
        dat = f.read()
    sequence = dat.split(">")[1:]
    for line in sequence:
        header,seq = line.split("\n",1)
        seq = seq.replace("\n","")
        FaD[header] = seq
    return(FaD)

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
def read_parameter_p3(filename):
    with open(filename, 'r') as f:
        lines = f.read().split('\n')
    parameters = {}
    for line in lines:
        if not line.startswith('#') and '=' in line:
            key, value = line.split('=')
            if '.' in value:
                parameters[key.strip()] = float(value)
            else:
                parameters[key.strip()] = int(value)
    return parameters

def CoreAndCloud(TotalGene,refID):
    dat = pd.DataFrame(TotalGene.apply(lambda x: [sum(x!=0)/TotalGene.shape[1], 
                                   np.mean(x), 
                                   np.std(x), 
                                   max(x)-min(x), 
                                   max(x)/(min(x)+1 if min(x)==0 else min(x)),
                                   sum(x)],axis=1))
    dat = dat[0].apply(pd.Series).reset_index()
    Count = pd.merge(dat,refID[["OGID","length"]],left_on="Orthogroup",right_on="OGID").drop(["OGID"],axis=1)
    Count.columns = ["Orthogroup","Perc","Mean","Var","Range","Fold","Total","length"]
    Count = Count.sort_values(["Perc","Mean","Var","Range","length"],ascending=[False,True,True,True,True])
    return(Count)


def Unique(TotalGene, refID):
    # Create a dictionary of unique columns for each row in TotalGene
    unique_cols = dict(TotalGene.apply(lambda x: x.to_numpy().nonzero()[0][0] if x.to_numpy().nonzero()[0].size == 1 else None, axis=1).dropna())
    unique_cols_dict = {x:TotalGene.columns.str.replace("_cds","")[y] for x,y in unique_cols.items()}


    # Create a DataFrame from the dictionary
    UniC = pd.DataFrame.from_dict(unique_cols_dict, orient='index', columns=['Cds'])
    UniC.index.name = 'Orthogroup'
    UniC.reset_index(inplace=True)

    # Merge with refID and drop unnecessary columns
    UniC = pd.merge(UniC, refID[['OGID', 'length']], left_on='Orthogroup', right_on='OGID').drop('OGID', axis=1)

    return(UniC)

def PrimerMake(Id,FaD,globe_args):
    pattern = re.compile('[^ATCG]')
    if not bool(pattern.search(FaD[Id])):
            seq_args = {
                    'SEQUENCE_ID':Id ,
                    'SEQUENCE_TEMPLATE': FaD[Id],
                    'SEQUENCE_INCLUDED_REGION': [0,len(FaD[Id])],
                    }
            print(Id,len(FaD[Id]))
            dat = primer3.bindings.designPrimers(seq_args,globe_args)
            df_tmp = primer_result(dat,Id)
    else:
            df_tmp = pd.DataFrame()
    return(df_tmp)


def primer_result(dat, ID):
    keys = list(dat.keys())
    left_L = [dat[k] for k in keys if k.startswith('PRIMER_LEFT_') and k.endswith('_SEQUENCE')]
    right_L = [dat[k] for k in keys if k.startswith('PRIMER_RIGHT_') and k.endswith('_SEQUENCE')]
    left_TM_L = [dat[k] for k in keys if k.startswith('PRIMER_LEFT_') and k.endswith('_TM')]
    right_TM_L = [dat[k] for k in keys if k.startswith('PRIMER_RIGHT_') and k.endswith('_TM')]
    left_GC_L = [dat[k] for k in keys if k.startswith('PRIMER_LEFT_') and k.endswith('_GC_PERCENT')]
    right_GC_L = [dat[k] for k in keys if k.startswith('PRIMER_RIGHT_') and k.endswith('_GC_PERCENT')]
    left_stable = [dat[k] for k in keys if k.startswith('PRIMER_LEFT_') and k.endswith('_END_STABILITY')]
    right_stable = [dat[k] for k in keys if k.startswith('PRIMER_RIGHT_') and k.endswith('_END_STABILITY')]
    product_size = [dat[k] for k in keys if k.startswith('PRIMER_PAIR_') and k.endswith('_PRODUCT_SIZE')]
    df_tmp = pd.DataFrame({
        'ID': [ID] * len(left_L),
        'LP': left_L,
        'RP': right_L,
        'TMl': left_TM_L,
        'TMr': right_TM_L,
        'GCl': left_GC_L,
        'GCr': right_GC_L,
        'Stablel': left_stable,
        'Stabler': right_stable,
        'ExceptSize': product_size,
    })
    return df_tmp

def PrimerDesign(TotalGene,refID,globe_args,Taxonomy,workPath):
    print("PrimerDesign Start!")
    pattern = re.compile('[^ATCG]')
    if Taxonomy == "Strain": 
        primerC = Unique(TotalGene,refID)
    else:
        CountC = CoreAndCloud(TotalGene,refID)
        primerC = CountC[CountC["Perc"] > 0.99]
    maxL = 1000
    minL = 150
    primerC = primerC.query(f"length > {minL} & length < {maxL}")
    FaD = Fa2D(os.path.join(workPath,"Total.OGID.fa"))
    TotalNum = len(primerC["Orthogroup"].drop_duplicates())

    logging.info(f"Total Primer Pre Orthogroups is {TotalNum}")
    logging.info("Primer Design Start!")

    results = Parallel(n_jobs=10)(delayed(PrimerMake)(Id,FaD,globe_args) for Id in primerC["Orthogroup"].drop_duplicates().tolist())
    # result = []
    # for Id in primerC["Orthogroup"].drop_duplicates().tolist():
    #     temp = PrimerMake(Id,FaD,globe_args)
    #     result.append(temp)
    df_primer = pd.concat(results).reset_index()
    df_primer.insert(2, "primer", df_primer.apply(lambda x: f"{x['ID']}_{x['index']}", axis=1))
    df_primer["Stable"] = df_primer[["Stablel","Stabler"]].mean(axis=1)
    df_primer = df_primer.drop(["Stabler","Stablel","index"], axis=1)


    lp_counts = df_primer["LP"].value_counts().to_dict()
    rp_counts = df_primer["RP"].value_counts().to_dict()
    df_primer["LP_count"] = df_primer["LP"].map(lp_counts)
    df_primer["RP_count"] = df_primer["RP"].map(rp_counts)
    logging.info("Primer Design complete!{} primers have been wriitern".format(df_primer.shape[0]))
    return(df_primer)

def SaveR(primer,bac,outpath):
    PrimerPath = os.path.join(outpath,"Working_Direct")
    os.makedirs(PrimerPath, exist_ok=True)

    with open(os.path.join(PrimerPath, f"{bac}_1.fa"), "w+") as leftP, \
        open(os.path.join(PrimerPath, f"{bac}_2.fa"), "w+") as rightP:

        for index, row in primer.iterrows():
            leftP.write(f">{row['primer']}_L\n{row['LP']}\n")
            rightP.write(f">{row['primer']}_R\n{row['RP']}\n")

    primer.drop(columns=['LP_count', 'RP_count'], inplace=True)
    primer.to_csv(os.path.join(PrimerPath, f"{bac}_primer_metadata.txt"), sep="\t", index=False)

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
    refseq = default["refseq"]
    globe_args = read_parameter_p3(default["p3_config"])
    Taxonomy = config["Taxonomy"]
    bac = config["Bacteria"].replace(" ","_")
    df_Meta = pd.read_table(os.path.join(config["Result"],bac,"%s_metadata.txt"%bac))
    for index,value in df_Meta.iterrows():
        workPath = os.path.join(config["Temp_file"],"%s/%s/OrthoFinder/Results_temp/Orthogroups"%(value["Genus"],value["Species"].replace(" ","_")))
        outPath = os.path.join(config["Result"],bac)

        TotalGene = pd.read_table(os.path.join(workPath,"Total.GeneCount.tsv"),index_col=0)
        refID = pd.read_table(os.path.join(workPath,"Total_sequences_refID2.txt"))

        Primerbac = PrimerDesign(TotalGene,refID,globe_args,Taxonomy,workPath)
        SaveR(Primerbac,bac,outPath)

if __name__ == "__main__":
    main()
