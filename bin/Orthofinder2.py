import pandas as pd
import os
import glob
import subprocess
import logging
import argparse
import re
import sys
from ete3 import NCBITaxa

# 设置日志记录器
logging.basicConfig(filename='orthofinder.log', level=logging.INFO, 
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


def parse_fasta_files(files):
    faD = {}
    faL = {}
    for file in files:
        P_dict, L_dict = Fa2DL(file)
        faD.update(P_dict)
        faL.update(L_dict)
    return(faD, faL)


def run_orthofinder(input_dir, use_docker=True, docker_image="davidemms/orthofinder:2.5.4"):
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-C", "--config",default=None,type=str,required=True,help="The input yaml file.")
    parser.add_argument("-B", "--bac",default=None,type=str,required=True,help="Input Bacteria Name.")
    parser.add_argument("-T", "--tax",default="Species",type=str,required=False,help="The Taxnomy level you want!defalut:Species.")
    parser.add_argument("-F", "--fasta",default=None,type=str,required=False,help="The FASTA file location.")
    parser.add_argument("-O", "--out",default=None,type=str,required=True,help="The oputput path.")
    parser.add_argument("-b", "--backend",default=None,type=str,required=False,help="The Previous Orthofinder result directory.")
    args = parser.parse_args()
    conF = args.config
    Taxonomy =  args.tax
    Path = args.out
    default = configP(conF)
    bacteria = args.bac
    backFolder = args.backend
    fastaD = args.fasta

    default = configP(conF)
    ncbi = NCBITaxa()
    #ncbi.update_taxonomy_database(default["taxonomy"])   
    Intens = pd.read_table(default["intense"])
    refseq = pd.read_table(default["refseq"],header=1)
    
    Temp = os.path.join(Path,"Temp")
    datMeta = IDtrans(ncbi,bacteria,Taxonomy,refseq,Intens)[2]
    for index,value in datMeta.iterrows():
        if not backFolder:
            newGe = value["Genus"]
            newSp = value["Species"].replace(" ","_").replace(".","")
            newPath = os.path.join(Temp,f"{newGe}/{newSp}")
            os.makedirs(newPath,exist_ok=True)
            run_orthofinder(newPath)
            resultPath = f"{newPath}/OrthoFinder/Results_temp/Orthogroups"
            PanG = PanGenome(resultPath,newPath)
        else:
            PanG = PanGenome(backFolder,fastaD)
        PanG.SaveR()

if __name__ =="__main__":
    main()