import os
import logging
import requests
import shutil
from ete3 import NCBITaxa
import pandas as pd
import numpy as np
import sys
import gzip
from joblib import Parallel, delayed
import argparse
import threading
import subprocess
import re
import primer3
from requests.exceptions import RequestException

lock = threading.Lock()

# 设置日志记录器
logging.basicConfig(filename='primerBac.log', level=logging.INFO, 
                    format='%(asctime)s %(levelname)s: %(message)s')

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


# ID转换函数
def TaxIDtrans(ncbi,Id=None,name=None):
    if name:
        Trans = ncbi.get_name_translator([name])
        tax = Trans[name][0] if Trans else None
    elif Id:
        Trans = ncbi.get_taxid_translator([Id])
        tax = Trans[Id] if Trans else None
    else:
        tax = None
    return(tax)


# 定义下载函数
def download_file(download_dir, gcf_url, file_name, base_url, retries=3):
    logging.info('开始下载文件: {}'.format(file_name))
    Id = gcf_url.split("/")[-1]
    url = base_url.format(gcf_url, Id)
    file_path = os.path.join(download_dir, file_name+'_cds.faa.gz')
    failed_dir =os.path.join(download_dir,"Failed_downloads")
    
    if not failed_dir:
        os.makedirs(failed_dir)
    
    for i in range(retries):
        try:
            # 下载文件
            response = requests.get(url, stream=True)
            with open(file_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=1024):
                    f.write(chunk)
            # 解压文件
            extracted_file_path = file_path[:-3]
            with gzip.open(file_path, 'rb') as f_in:
                with open(extracted_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            # 删除原文件
            os.remove(file_path)
            logging.info(f'{file_name}下载完成')
            return True
        except RequestException as e:
            logging.warning(f'{file_name}下载失败，错误信息：{str(e)}')
            if i == retries-1:
                failed_file_path = os.path.join(failed_dir, file_name)
                os.rename(file_path, failed_file_path)
                logging.warning(f'{file_name}已添加到失败下载目录')
                return False
            
def IDtrans(ncbi,bac,tax,refseq,intens):
    taxid = TaxIDtrans(ncbi,name=bac)
    logging.info('ID Trans Complete: {} to {}'.format(bac,taxid))
    if taxid:
        rankD = dict([(y.capitalize(),x) for x,y in ncbi.get_rank(ncbi.get_lineage(taxid)).items()])
        if tax == "Genus":
            Genus = TaxIDtrans(ncbi,Id=rankD["Genus"])
            SpeL = intens.loc[intens["Genus"] == Genus,"Species"].drop_duplicates()
            taxidL = list(map(lambda x:TaxIDtrans(ncbi,name=x),SpeL))
            datDown = refseq[refseq["species_taxid"].isin(taxidL)]
            datMeta = pd.DataFrame([taxidL,[Genus]*len(taxidL),list(SpeL),list(map(lambda x:datDown[datDown["species_taxid"]==x].shape[0],taxidL))],
                                index=["TaxID","Genus","Species","RefSeqNum"]).T
        else:
            datDown = refseq[refseq["species_taxid"] == taxid]
            datMeta = pd.DataFrame([taxid,TaxIDtrans(ncbi,Id=rankD["Genus"]),TaxIDtrans(ncbi,Id=rankD["Species"]),len(datDown)],index=["TaxID","Genus","Species","RefSeqNum"]).T
    else:
        logging.info("The bacteriad name  %s is not right!"%bac)
        sys.exit()
    return(rankD,datDown,datMeta)

            
def DownFile(rankD,datDown,datMeta,ncbi,bac,tax,temp,thread):
    Genus = TaxIDtrans(ncbi,Id=rankD["Genus"])
    base_url = "{}/{}_cds_from_genomic.fna.gz"
    if len(datDown["ftp_path"]) > 500:
        print("The bacteriad  %s sample is too large!"%bac)
        sys.exit()
    outPath  = []
    if tax == "Genus":
        for i,v in datMeta.iterrows():
            dat_temp = datDown.query("species_taxid=={}".format(v["TaxID"]))
            gcf_urls = dat_temp["ftp_path"].tolist()
            file_names = [x.split("/")[-1].split(".")[0] for x in gcf_urls]
            newPath = "%s/%s/%s"%(temp,Genus,v["Species"].replace(" ","_"))
            outPath.append(newPath)
            if not os.path.exists(newPath):
                os.makedirs(newPath)
            result = Parallel(n_jobs=thread)(delayed(download_file)(newPath,gcf_url,file_name,base_url) for gcf_url, file_name in zip(gcf_urls,file_names))
                # 输出下载结果                 
            logging.info('Species {} files is all completed!'.format(str(v["Species"]))) if all(result) else  print('存在下载失败的文件，请查看日志文件')
    else:
        gcf_urls = datDown["ftp_path"].tolist()
        file_names = [x.split("/")[-1].split(".")[0] for x in gcf_urls]
        newPath = "%s/%s/%s"%(temp,Genus,bac.replace(" ","_"))
        outPath.append(newPath)
        if not os.path.exists(newPath):
                os.makedirs(newPath)
        result = Parallel(n_jobs=thread)(delayed(download_file)(newPath,gcf_url,file_name,base_url) for gcf_url,file_name in zip(gcf_urls,file_names))
        logging.info('Species {} files is all completed!'.format(datMeta["Species"].str.cat())) if all(result) else  print('存在下载失败的文件，请查看日志文件')

    return(outPath)

def Fa2DL(InF):
    P_dic = {}
    with open(InF,"r") as f:
        dat = f.readlines()
    for line in dat:
        line = line.strip()
        if line.startswith(">"):
            nam = line.strip(">").split(" ")[0]
            P_dic[nam] = ""
        else:
            P_dic[nam] += line
    L_dic = {k: len(v) for k, v in P_dic.items()}
    return(P_dic, L_dic)

def parse_fasta_files(files):
    faD = {}
    faL = {}
    for file in files:
        P_dict, L_dict = Fa2DL(file)
        faD.update(P_dict)
        faL.update(L_dict)
    return(faD, faL)


def run_orthofinder(input_dir,use_docker=True, docker_image="davidemms/orthofinder:2.5.4"):
    logger = logging.getLogger(__name__)

    cmd = ["orthofinder", "-f", input_dir,"-n","temp","-S","diamond","-M","msa","-os", "-t","208","-a","208","-d"]
    if use_docker:
        import os
        uid = os.getuid()
        gid = os.getgid()
        cmd = ["docker", "run", "--ulimit","nofile=1000000:1000000","--rm","-u", f"{uid}:{gid}", "-v", f"{input_dir}:{input_dir}:Z", docker_image] + cmd

    try:
        subprocess.run(cmd, check=True)
        logger.info("Orthofinder finished running successfully")
    except subprocess.CalledProcessError as e:
        logger.error(f"Orthofinder failed with error: {e}")
        raise e
    
class PanGenome(object):
    def __init__(self,workpath,fastaD):
        self.workpath = workpath
        self.fastaD = fastaD
        self.datAll = self.GeneCountTable()
        self.datID = self.refID()
        self.FaD = self.get_sequence_data()[0]
        self.length = self.get_sequence_data()[1]

    
    def GeneCountTable(self):
        df_Genecount = pd.read_table(os.path.join(self.workpath, "Orthogroups.GeneCount.tsv"))
        df_UnassignedGenes = pd.read_table(os.path.join(self.workpath, "Orthogroups_UnassignedGenes.tsv"), index_col=0)

        df_Genecount.columns = df_Genecount.columns.str.replace("_protein|_cds", "")
        df_UnassignedGenes.columns = df_UnassignedGenes.columns.str.replace("_protein|_cds", "")
        df_UnassignedGenes = 1 - df_UnassignedGenes.isna().astype(int)

        df_UnassignedGenes = df_UnassignedGenes.reset_index()
        df_Genecount = df_Genecount.iloc[:, :-1]

        df_all = pd.concat([df_Genecount, df_UnassignedGenes])
        return(df_all)
    
    def refID(self):
        with open(os.path.join(self.workpath, "Orthogroups.txt"), "r") as f:
            dat_Orth = f.readlines()
        with open(os.path.join(self.workpath, "Orthogroups_UnassignedGenes.tsv"), "r") as f:
            dat_UnassignedGenes = f.readlines()

        refIDL = [[tmp1[0], n] for i in dat_Orth
                                for tmp1 in [i.strip().split(":")]
                                for n in tmp1[1].split(" ") if n]
        refIDL.extend([[elem.strip() for elem in sublist if elem.strip()] 
                        for sublist in [x.strip().split("\t") for x in dat_UnassignedGenes[1:]]])

        datID = pd.DataFrame(refIDL, columns=["OGID", "GeneID"]).drop_duplicates()
        return(datID)
    

    def get_sequence_data(self):
        pattern = re.compile(r"\.(faa|fasta|pep|fa)$")
        files = [os.path.join(self.fastaD,f) for f in os.listdir(self.fastaD) if pattern.search(f)]
        FaD,FaL = parse_fasta_files(files)
        datLength = pd.DataFrame(FaL.items(), columns=["GeneID", "length"])
        return(FaD, datLength)

    def SaveR(self):
        # Write dataframes to files
        self.datAll.to_csv(os.path.join(self.workpath, "Total.GeneCount.tsv"), sep="\t", index=False)
        self.datID.to_csv(os.path.join(self.workpath, "Total_sequences_refID2oriID.txt"), sep="\t", index=False)
        self.length.to_csv(os.path.join(self.workpath, "PepSequencesLength.txt"), sep="\t", index=False)

        # Merge dataframes and sort by OGID and length
        dat_full = pd.merge(self.datID, self.length).sort_values(["OGID", "length"], ascending=(True, False))

        # Write merged dataframe to file
        dat_full.to_csv(os.path.join(self.workpath, "Total_sequences_refID2oriID2.txt"), sep="\t", index=False)

        # Get reference sequence for each OGID
        data_ref = dat_full.loc[dat_full.groupby("OGID")["length"].idxmin()]
        dat_ref_anno = dat_full.loc[dat_full.groupby("OGID")["length"].idxmax()]
        dat_ref_anno.to_csv(os.path.join(self.workpath,"Total_sequences_refID.txt"),sep="\t",index=False)
        data_ref.to_csv(os.path.join(self.workpath,"Total_sequences_refID2.txt"),sep="\t",index=False)

        # Write reference sequence IDs to file
        data_ref["GeneID"].to_csv(os.path.join(self.workpath, "GeneID.list"), sep="\t", index=False, header=None)

        # Write reference sequences to files
        with open(os.path.join(self.workpath, "Total.oriID.fa"), "w+") as out:
            for i in data_ref["GeneID"].to_list():
                out.writelines("%s\t%s\n" % (i, self.FaD[i]))

        with open(os.path.join(self.workpath, "Total.OGID.fa"), "w+") as out:
            for index, row in data_ref.iterrows():
                out.writelines(">%s\n%s\n" % (row["OGID"], self.FaD[row["GeneID"]]))

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
        'Orthogroup': [ID] * len(left_L),
        'LP': left_L,
        'RP': right_L,
        'TMl': left_TM_L,
        'TMr': right_TM_L,
        'GCl': left_GC_L,
        'GCr': right_GC_L,
        'Stable': np.mean(left_stable+right_stable),
        'ExceptSize': product_size,
    })
    return df_tmp

def RevRatioAndGC(seq1, seq2):
    # 计算反向互补序列
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    revcomp1 = ''.join([complement[base] for base in seq1[::-1]])
    revcomp2 = ''.join([complement[base] for base in seq2[::-1]])
    
    # 计算反向互补比例
    revcomp_ratio1 = sum([1 for base1, base2 in zip(seq1, revcomp1) if base1 == base2]) / len(seq1)
    revcomp_ratio2 = sum([1 for base1, base2 in zip(seq2, revcomp2) if base1 == base2]) / len(seq2)
    revcomp_ratio = max(revcomp_ratio1,revcomp_ratio2)*100

    ## 计算引物间互补比例
    comp2 = ''.join([complement[base] for base in seq2])
    comp_count = 0
    max_comp_count = 0
    for base1, base2 in zip(seq1, comp2):
        if base1 == base2:
            comp_count +=1
        else:
            comp_count =0
        max_comp_count = max(comp_count,max_comp_count)
        
    # 计算最大连续GC数量
    gc_count = 0
    max_gc_count = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 == 'G' or base1 == 'C':
            gc_count += 1
        else:
            gc_count = 0
        max_gc_count = max(max_gc_count, gc_count)
        
        if base2 == 'G' or base2 == 'C':
            gc_count += 1
        else:
            gc_count = 0
        max_gc_count = max(max_gc_count, gc_count)
    
    return revcomp_ratio,max_comp_count, max_gc_count


def Bestprimer(dat,default):

    Continue_GC_Num_Max = default["Continue_GC_Num_Max"]
    Continue_GC_Num_Min = default["Continue_GC_Num_Min"]
    Reverse_Ration_Max = default["Reverse_Ration_Max"]
    Reverse_Ration_Min = default["Reverse_Ration_Min"]
    Primer_Reverse_Num = default["Primer_Reverse_Num"]
    Primer_size_Min = default["Primer_size_Min"]
    Primer_size_Max = default["Primer_size_Max"]

    rev_ratio_and_gc  = dat.apply(lambda row:RevRatioAndGC(row["LP"],row["RP"]),axis=1)
    dat["RevRatio"] = [item[0] for item in rev_ratio_and_gc]
    dat["RevComp"] = [item[1] for item in rev_ratio_and_gc]
    dat["GCNum"] = [item[2] for item in rev_ratio_and_gc]
    dat.insert(2, "primer", dat.apply(lambda x: f"{x['Orthogroup']}_{x['index']}", axis=1))
    dat = dat.drop(["index"], axis=1)
    lp_counts = dat["LP"].value_counts().to_dict()
    rp_counts = dat["RP"].value_counts().to_dict()
    dat["LP_count"] = dat["LP"].map(lp_counts)
    dat["RP_count"] = dat["RP"].map(rp_counts)
    best_primer = dat.query(f"RevRatio >= {Reverse_Ration_Min} & RevRatio < {Reverse_Ration_Max} &\
                    GCNum >={Continue_GC_Num_Min} & GCNum < {Continue_GC_Num_Max} &\
                    ExceptSize >= {Primer_size_Min} & ExceptSize < {Primer_size_Max} &\
                    RevComp < {Primer_Reverse_Num}")
    return(best_primer)

def PrimerDesign(TotalGene,refID,globe_args,Taxonomy,workPath,default,threads):
    CountC = CoreAndCloud(TotalGene,refID)
    UniC = Unique(TotalGene,refID)
    CountC["Type"] = CountC["Perc"].apply(lambda x:"Core" if x>0.95 else None)
    CountC.loc[CountC["Orthogroup"].isin(UniC["Orthogroup"]),"Type"] = "Unique"
    CountC["Type"].fillna("Accessory",inplace=True)
    primerC = CountC.query('Type == "Unique"') if Taxonomy == "Strain" else CountC.query('Type == "Core"')

    maxL = 1000
    minL = 150
    primerC = primerC.query(f"length > {minL} & length < {maxL}")
    FaD = Fa2D(os.path.join(workPath,"Total.OGID.fa"))

    results = Parallel(n_jobs=threads)(delayed(PrimerMake)(Id,FaD,globe_args) for Id in primerC["Orthogroup"].drop_duplicates().tolist())
    # result = []
    # for Id in primerC["Orthogroup"].drop_duplicates().tolist():
    #     temp = PrimerMake(Id,FaD,globe_args)
    #     result.append(temp)
    df_primer = pd.concat(results).reset_index()
    best_primer = Bestprimer(df_primer,default)
    return(CountC,best_primer)

def SaveR(primer,bac,PrimerPath):
    os.makedirs(PrimerPath, exist_ok=True)

    with open(os.path.join(PrimerPath, f"{bac}_1.fa"), "w+") as leftP, \
        open(os.path.join(PrimerPath, f"{bac}_2.fa"), "w+") as rightP:

        for index, row in primer.iterrows():
            leftP.write(f">{row['primer']}_L\n{row['LP']}\n")
            rightP.write(f">{row['primer']}_R\n{row['RP']}\n")

    primer.drop(columns=['LP_count', 'RP_count'], inplace=True)
    primer.to_csv(os.path.join(PrimerPath, f"{bac}_primer_metadata.txt"), sep="\t", index=False)

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

def PrimerBac(value,default,globe_args,Taxonomy,bac,bowtieIn,id2sp,BacNum,threads):
    newG = value["Genus"]
    newSp = value["Species"].replace(" ","_").replace(".","")
    workPath = os.path.join(default["Temp_file"],f"{newG}/{newSp}/OrthoFinder/Results_temp/Orthogroups")
    outPath = os.path.join(default["Result"],f"{newG}/{newSp}")
    PrimerPath = os.path.join(outPath,"Working_Direct")

    TotalGene = pd.read_table(os.path.join(workPath,"Total.GeneCount.tsv"),index_col=0)
    refID = pd.read_table(os.path.join(workPath,"Total_sequences_refID2.txt"))
    
    logging.info("{} Primer Design Started!".format(value["Species"]))
    PrimerC,Primerbac = PrimerDesign(TotalGene,refID,globe_args,Taxonomy,workPath,default,threads)
    logging.info("{} Primer Design complete!{} primers have been wriitern".format(value["Species"],Primerbac.shape[0]))
    SaveR(Primerbac,bac,PrimerPath)

    bt = Bowtie(bowtieIn,PrimerPath,outPath,bac)
    bt.run()
    bacteriaName = value["Genus"] if default["Taxonomy"]=="Genus" else value["Species"]
    dat = SpeciesR(os.path.join(outPath,f"{bac}.sam"),id2sp,bacteriaName,BacNum)
    dat = dat.merge(Primerbac,on="primer",how="left")
    Primerbac = dat.loc[dat.groupby("Orthogroup")["Stable"].idxmax(),:]
    PrimerC.to_csv(os.path.join(outPath,f"{newSp}_pangenome.txt"),index=False,sep="\t")
    return(Primerbac)




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
    parser.add_argument(
        "-b",
        "--bacteria",
        default=None,
        type=str,
        required=True,
        help="The Bacteria Name you want.",
    )
    args = parser.parse_args()
    conF = args.config
    bacteria = args.bacteria
    default = configP(conF)
    ncbi = NCBITaxa()
    #ncbi.update_taxonomy_database(default["taxonomy"])   
    Intens = pd.read_table(default["intense"])
    refseq = pd.read_table(default["refseq"],header=1)
    
    Temp = default["Temp_file"]
    threads = default["Thread"]
    globe_args = read_parameter_p3(default["p3_config"])
    Taxonomy = default["Taxonomy"]
    df_intense = pd.read_table(default["Taxanno"])
    id2sp = dict(zip(df_intense["Gene_ID"],df_intense["Species"]))
    Ge2Spnum =df_intense[["Genus","Species"]].drop_duplicates()["Genus"].value_counts().to_dict()
    Sp2StNum = df_intense["Species"].value_counts().to_dict()
    BacNum = Ge2Spnum if default["Taxonomy"] == "Genus" else Sp2StNum
    bowtieIn = default["BowIndex"]


    rankD,datDown,datMeta = IDtrans(ncbi,bacteria,Taxonomy,refseq,Intens)
    datMeta = datMeta.query("RefSeqNum > 1")
    genus = TaxIDtrans(ncbi,Id=rankD["Genus"])
    species = TaxIDtrans(ncbi,Id=rankD["Species"]).replace(" ","_").replace(".","")
    OutBac =  genus if default["Taxonomy"] == "Genus" else species.strip().replace(" ","_").replace(".","")
    ResultPath = os.path.join(default["Result"],f"{genus}/{species}")
    if not os.path.exists(ResultPath):
        os.makedirs(ResultPath)
    outPath = os.path.join(default["Result"],f"{genus}") if default["Taxonomy"] == "Genus" else os.path.join(default["Result"],f"{genus}/{species}")
    datDown.to_csv(os.path.join(outPath,"%s_download.txt"%OutBac),sep="\t",index=False)
    datMeta.to_csv(os.path.join(outPath,"%s_metadata.txt"%OutBac),sep="\t",index=False)
    
    FilePath  = DownFile(rankD,datDown,datMeta,ncbi,bacteria,Taxonomy,Temp,threads)
    for path in FilePath:
        run_orthofinder(path)
        resultPath = "%s/OrthoFinder/Results_temp"%path
        WorkDirectory = os.path.join(resultPath,"Orthogroups")
        PanG = PanGenome(WorkDirectory,path)
        PanG.SaveR()

    results = [PrimerBac(value,default,globe_args,Taxonomy,species,bowtieIn,
                                         id2sp,BacNum,threads) for index,value in datMeta.iterrows()]
    datPrimer = pd.concat(results)
    FinalP = os.path.join(outPath,f"{genus}.primer.txt") if default["Taxonomy"] == "Genus" else os.path.join(outPath,f"{species}.primer.txt")
    datPrimer.sort_values("Stable",ascending = False).to_csv(FinalP,sep="\t",index=False)


if __name__ == "__main__":
    main()