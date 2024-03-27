import os
import subprocess
import re
from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse

def create_host_database(pathogen_files, host_files, hostdb_file):
    # Concatenate the pathogen and host proteomes into separate files
    cat_pathogen_command = "cat " + " ".join(pathogen_files) + " > pathogens.faa"
    cat_host_command = "cat " + " ".join(host_files) + " > host.faa"
    subprocess.call(cat_pathogen_command, shell=True)
    subprocess.call(cat_host_command, shell=True)

    # Create a Diamond database from the merged host proteome file
    if not os.path.isfile(hostdb_file + ".dmnd"):
        makedb_command = f"diamond makedb --in host.faa --db {hostdb_file}"
        subprocess.call(makedb_command, shell=True)
    else:
        print(f"Host database file {hostdb_file}.dmnd already exists.")

def fix_fasta_file(filename):
    fixed_filename = f"{filename[:-5]}_fixed.fasta"
    subprocess.run(["seqret", "-sequence", filename, "-osformat", "fasta", "-auto", "-outseq", fixed_filename])
    return fixed_filename

def merge_domains(domains, threshold=30):
    """
    Merge overlapping or close domains based on a distance threshold
    """
    merged_domains = []
    domains = sorted(domains, key=lambda x: x[1])
    current_domain = domains[0]
    for domain in domains[1:]:
        # If the domains overlap or are within the distance threshold, merge them
        if domain[1] < current_domain[2] or domain[1] - current_domain[2] <= threshold:
            current_domain = (current_domain[0], current_domain[1], max(domain[2], current_domain[2]))
        else:
            merged_domains.append(current_domain)
            current_domain = domain
    merged_domains.append(current_domain)
    return merged_domains

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pathogen Protein Analysis Pipeline")

    # Step 1 arguments
    parser.add_argument("--host_proteomes", nargs="+", help="Path to the host proteome file", required=True)
    parser.add_argument("--pathogen_proteomes", nargs="+", help="List of pathogen proteome files", required=True)
    parser.add_argument("--step1_output", default="step1.nonhost.faa", help="Output file for non-host pathogen protein sequences")
    parser.add_argument("--hostdb_file", default="hostdb", help="Name of the host database file")

    # Step 2 arguments
    parser.add_argument("--deg_database", default="DEG20.faa", help="Path to the DEG database file", required=True)
    parser.add_argument("--step2_output", default="step2.nohost.deg.faa", help="Output file for DEG-similar pathogen protein sequences")

    # Step 3 arguments
    parser.add_argument("--go_annotations", default="go-basic.obo", help="Path to the gene ontology annotations file", required=True)
    parser.add_argument("--step3_output", default="step3.nohost.deg.go.faa", help="Output file for GO annotated pathogen protein sequences")

    # Step 4 arguments
    parser.add_argument("--pfam_database", default="Pfam-A.hmm", help="Path to the Pfam database directory", required=True)
    parser.add_argument("--step4_output", default="step4.domain.faa", help="Output file for domain annotated pathogen protein sequences")

    # Step 5 arguments
    parser.add_argument("--step5_domain_output", default="step5.nohost.domain.faa", help="Output file for non-host domain pathogen protein sequences")
    parser.add_argument("--step5_fullpro_output", default="step5.nohost.deg.go.domain.faa", help="Output file for non-host domain pathogen full protein sequences")

    # Step 6 arguments
    parser.add_argument("--step6_output", default="step6.nohost.deg.go.domain.drugable.faa", help="Output file for drugable pathogen protein sequences")
    parser.add_argument("--step6_csvoutput", default="final.result.csv", help="Output file for total analysis result in csv file")

    args = parser.parse_args()

    pathogens_file = "pathogens.faa"
    deg_file = args.deg_database
    pfam_db = args.pfam_database
    goann_file = args.go_annotations

    step1_file = args.step1_output
    step2_file = args.step2_output
    step3_file = args.step3_output
    step4_file = args.step4_output
    step5_file_1 =  args.step5_domain_output
    step5_file_2 = args.step5_fullpro_output
    step6_file = args.step6_output
    final_file = args.step6_csvoutput

    create_host_database(args.pathogen_proteomes, args.host_proteomes, args.hostdb_file)

    step0_data = []
    with open(pathogens_file) as input_file:  
        # Process each line in the input file
        for line in input_file:
            if line.startswith(">"):
                # Extract the sequence ID
                seq_id, protein_name = line[1:].split(None, 1)
                protein_name = protein_name.split("OS=")[0].strip()
                protein_name = protein_name.strip('\"')
                step0_data.append({"ID": seq_id, "Protein Name": protein_name})

    # Create a DataFrame from the data list
    df_protein = pd.DataFrame(step0_data)

    ################## STEP1 Identify non-host pathogenic protein sequences ############################

    # Perform a Diamond BLASTp search using the merged host proteome database
    blastp_output_file = "inhost.txt"
    if os.path.isfile(blastp_output_file):
        os.remove(blastp_output_file)

    subprocess.call(f"diamond blastp -q {pathogens_file} -d hostdb --outfmt 6 qseqid sseqid pident evalue --id 90 --subject-cover 80 --query-cover 70 --max-target-seqs 1 --out {blastp_output_file}", shell=True)

    # Parse the BLASTp output to extract the inhost IDs
    inhost_ids = set()
    step1_data = []
    if os.path.isfile(blastp_output_file):
        with open(blastp_output_file) as f:
            for line in f:
                inhost_id = line.split("\t")[0].split()[0]
                host_sim = line.split("\t")[2].split()[0]
                inhost_ids.add(inhost_id)
                step1_data.append({inhost_id: host_sim})
    else:
        print(f"No Diamond BLASTp output file found: {blastp_output_file}")

    with open(pathogens_file) as f:
        with open(step1_file, "w") as out_file:
            seq_id = None
            seq = ""
            for line in f:
                if line.startswith(">"):
                    if seq_id is not None and seq_id.split()[0] not in inhost_ids:
                        out_file.write(f">{seq_id}\n{seq}\n")
                    seq_id = line.strip()[1:]
                    seq = ""
                else:
                    seq += line.strip()
            if seq_id is not None and seq_id.split()[0] not in inhost_ids:
                out_file.write(f">{seq_id}\n{seq}\n")

    step1_data = [(k, float(v)) for d in step1_data for k, v in d.items()]
    df_step1 = pd.DataFrame.from_dict(step1_data)
    df_step1.columns = ["ID", "host.sim"]

    df_total = pd.merge(df_protein, df_step1, on="ID", how="left")
    df_total["host.sim"].fillna(0, inplace=True)

    ################## STEP2 Identify pathogenic protein sequences similar to DEG ############################

    deg_file_fixed = fix_fasta_file(deg_file)

    # Perform a Diamond BLASTp search using the merged DEG proteome database
    degdb_file = "degdb"
    if not os.path.isfile(degdb_file + ".dmnd"):
        subprocess.call(f"diamond makedb --in {deg_file_fixed} --db {degdb_file}", shell=True)
    else:
        print(f"DEG database file {degdb_file}.dmnd already exists.")

    blastp_output_file = "indeg.txt"
    if os.path.isfile(blastp_output_file):
        os.remove(blastp_output_file)

    subprocess.call(f"diamond blastp -q {step1_file} -d degdb --outfmt 6 qseqid sseqid pident evalue --id 50 --subject-cover 80 --max-target-seqs 1 --out {blastp_output_file}", shell=True)

    # Parse the BLASTp output to extract the indeg IDs
    indeg_ids = set()
    step2_data = []
    if os.path.isfile(blastp_output_file):
        with open(blastp_output_file) as f:
            for line in f:
                indeg_id = line.split("\t")[0].split()[0]
                deg_sim = line.split("\t")[2].split()[0]
                indeg_ids.add(indeg_id)
                step2_data.append({indeg_id: deg_sim})
    else:
        print(f"No Diamond BLASTp output file found: {blastp_output_file}")

    with open(pathogens_file) as seq_file:
        seq_id = None
        seq = ""
        with open(step2_file, "w") as out_file:
            for line in seq_file:
                if line.startswith(">"):
                    if seq_id is not None and seq_id.split()[0] in indeg_ids:
                        out_file.write(f">{seq_id}\n{seq}\n")
                    seq_id = line.strip()[1:]
                    seq = ""
                else:
                    seq += line.strip()
            if seq_id is not None and seq_id.split()[0] in indeg_ids:
                out_file.write(f">{seq_id}\n{seq}\n")

    step2_data = [(k, float(v)) for d in step2_data for k, v in d.items()]
    df_step2 = pd.DataFrame.from_dict(step2_data)
    df_step2.columns = ["ID", "deg.sim"]

    df_total = pd.merge(df_total, df_step2, on="ID", how="left")
    df_total.loc[df_total["host.sim"] > 0, "deg.sim"] = "uncheck"
    df_total["deg.sim"].fillna(0, inplace=True)

    ################## STEP3 Conduct gene ontology annotations on the sequences ############################

    target_go_terms = [
        {"id": "GO:0003824", "name": "catalytic activity", "namespace": "molecular_function"},
        {"id": "GO:0016740", "name": "transferase activity", "namespace": "molecular_function"},
        {"id": "GO:0005215", "name": "transporter activity", "namespace": "molecular_function"},
        {"id": "GO:0008152", "name": "metabolic process", "namespace": "biological_process"},
        {"id": "GO:0007047", "name": "cellular homeostasis", "namespace": "biological_process"},
        {"id": "GO:0055114", "name": "oxidation-reduction process", "namespace": "biological_process"},
        {"id": "GO:0009058", "name": "biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0009253", "name": "peptidoglycan biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0009252", "name": "peptidoglycan turnover", "namespace": "biological_process"},
        {"id": "GO:0009405", "name": "pathogenesis", "namespace": "biological_process"},
        {"id": "GO:0009406", "name": "toxin metabolic process", "namespace": "biological_process"},
        {"id": "GO:0016787", "name": "hydrolase activity", "namespace": "molecular_function"},
        {"id": "GO:0003677", "name": "DNA binding", "namespace": "molecular_function"},
        {"id": "GO:0003700", "name": "DNA-binding transcription factor activity", "namespace": "molecular_function"},
        {"id": "GO:0030246", "name": "carbohydrate binding", "namespace": "molecular_function"},
        {"id": "GO:0005576", "name": "extracellular region", "namespace": "cellular_component"},
        {"id": "GO:0005618", "name": "cell wall", "namespace": "cellular_component"},
        {"id": "GO:0009279", "name": "cell outer membrane", "namespace": "cellular_component"},
        {"id": "GO:0009404", "name": "toxin activity", "namespace": "molecular_function"},
        {"id": "GO:0008236", "name": "serine-type peptidase activity", "namespace": "molecular_function"},
        {"id": "GO:0030430", "name": "host cell cytoplasm", "namespace": "cellular_component"},
        {"id": "GO:0016020", "name": "membrane", "namespace": "cellular_component"},
        {"id": "GO:0006355", "name": "regulation of transcription, DNA-templated", "namespace": "biological_process"},
        {"id": "GO:0043066", "name": "negative regulation of apoptotic process", "namespace": "biological_process"},
        {"id": "GO:0051301", "name": "cell division", "namespace": "biological_process"},
        {"id": "GO:0016620", "name": "oxidoreductase activity", "namespace": "molecular_function"},
        {"id": "GO:0006508", "name": "proteolysis", "namespace": "biological_process"},
        {"id": "GO:0004672", "name": "protein kinase activity", "namespace": "molecular_function"},
        {"id": "GO:0006633", "name": "fatty acid biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0004152", "name": "dihydroorotate dehydrogenase activity", "namespace": "molecular_function"},
        {"id": "GO:0006164", "name": "purine nucleotide biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0008299", "name": "isoprenoid biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0005524", "name": "ATP binding", "namespace": "molecular_function"},
        {"id": "GO:0009055", "name": "electron transfer activity", "namespace": "molecular_function"},
        {"id": "GO:0006412", "name": "translation", "namespace": "biological_process"},
        {"id": "GO:0019856", "name": "pyrimidine nucleobase biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0033014", "name": "tetrapyrrole biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0003824", "name": "catalytic activity", "namespace": "molecular_function"},
        {"id": "GO:0004372", "name": "glycine hydroxymethyltransferase activity", "namespace": "molecular_function"},
        {"id": "GO:0016491", "name": "oxidoreductase activity", "namespace": "molecular_function"},
        {"id": "GO:0022857", "name": "transmembrane transporter activity", "namespace": "molecular_function"},
        {"id": "GO:0006096", "name": "glycolytic process", "namespace": "biological_process"},
        {"id": "GO:0016746", "name": "acyltransferase activity", "namespace": "molecular_function"},
        {"id": "GO:0008757", "name": "S-adenosylmethionine-dependent methyltransferase activity", "namespace": "molecular_function"},
        {"id": "GO:0004129", "name": "cytochrome-c oxidase activity", "namespace": "molecular_function"},
        {"id": "GO:0004422", "name": "hypoxanthine phosphoribosyltransferase activity", "namespace": "molecular_function"},
        {"id": "GO:0016114", "name": "terpenoid biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0005525", "name": "GTP binding", "namespace": "molecular_function"},
        {"id": "GO:0005975", "name": "carbohydrate metabolic process", "namespace": "biological_process"},
        {"id": "GO:0009239", "name": "enterobactin biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0008654", "name": "phospholipid biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0019752", "name": "carboxylic acid metabolic process", "namespace": "biological_process"},
        {"id": "GO:0071949", "name": "FAD binding", "namespace": "molecular_function"},
        {"id": "GO:0008610", "name": "lipid biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0007165", "name": "signal transduction", "namespace": "biological_process"},
        {"id": "GO:0006633", "name": "fatty acid biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0052636", "name": "arabinosyltransferase activity", "namespace": "molecular_function"},
        {"id": "GO:0003677", "name": "DNA binding", "namespace": "molecular_function"},
        {"id": "GO:0051920", "name": "peroxiredoxin activity", "namespace": "molecular_function"},
        {"id": "GO:0050661", "name": "NADP binding", "namespace": "molecular_function"},
        {"id": "GO:0006522", "name": "alanine metabolic process", "namespace": "biological_process"},
        {"id": "GO:0006221", "name": "pyrimidine nucleotide biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0016740", "name": "transferase activity", "namespace": "molecular_function"},
        {"id": "GO:0008830", "name": "dTDP-4-dehydrorhamnose 3,5-epimerase activity", "namespace": "molecular_function"},
        {"id": "GO:0008610", "name": "lipid biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0005506", "name": "iron ion binding", "namespace": "molecular_function"},
        {"id": "GO:0016747", "name": "acyltransferase activity, transferring groups other than amino-acyl groups", "namespace": "molecular_function"},
        {"id": "GO:0016757", "name": "glycosyltransferase activity", "namespace": "molecular_function"},
        {"id": "GO:0009435", "name": "NAD biosynthetic process", "namespace": "biological_process"},
        {"id": "GO:0019134", "name": "glucosamine-1-phosphate N-acetyltransferase activity", "namespace": "molecular_function"},
        {"id": "GO:0006799", "name": "polyphosphate biosynthetic process", "namespace": "biological_process"}
    ]

    subprocess.call(f"interproscan.sh -cpu 30 -i {step2_file} -f tsv -dp -goterms --iprlookup --pathways --appl Pfam -o {step3_file}.tsv", shell=True)

    go_terms = {}
    with open(goann_file, "r") as f:
        current_term = None
        for line in f:
            line = line.strip()
            if line == "[Term]":
                if current_term is not None:
                    go_terms[current_term["id"]] = current_term
                current_term = {}
            elif line.startswith("id: "):
                current_term["id"] = line[4:]
            elif line.startswith("name: "):
                current_term["name"] = line[6:]
            elif line.startswith("namespace: "):
                current_term["namespace"] = line[11:]
            elif line.startswith("is_a: "):
                if "is_a" not in current_term:
                    current_term["is_a"] = []
                current_term["is_a"].append(line[6:].split()[0])

    with open(step3_file, "w") as out:
        step3_data = []
        with open(f"{step3_file}.tsv", "r") as f:
            seq_dict = SeqIO.to_dict(SeqIO.parse(step2_file, "fasta"))
            merged_sequences = {}
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 14:
                    continue
                if fields[3] != "Pfam":
                    continue
                seq_id = fields[0]
                if seq_id not in seq_dict:
                    continue
                seq = seq_dict[seq_id]

                go_matches = re.findall(r"GO:(\d+)", fields[13])
                go_matches = ["GO:" + go_id for go_id in go_matches]            
                if not go_matches:
                    continue
                
                go_terms_desc = []
                go_ids = []
                for go_id in go_matches:
                    if go_id in go_terms:
                        go_term = go_terms[go_id]
                        for target_go_term in target_go_terms:
                            if go_term["id"] == target_go_term["id"] or go_term["name"] == target_go_term["name"]:
                                go_terms_desc.append(f"{go_term['name']} ({go_term['namespace']})")
                                go_ids.append(go_id)
                                break
                if not go_terms_desc:
                    continue
                
                go_terms_desc_str = ";".join(go_terms_desc)
                go_ids_str = ";".join(go_ids)

                if seq.id in merged_sequences:
                    merged_sequences[seq.id].description = go_terms_desc_str
                else:
                    merged_seq = seq_dict[seq.id]
                    merged_seq.description = go_terms_desc_str
                    merged_sequences[seq.id] = merged_seq

                step3_data.append({"ID": seq_id, "GO ID": go_ids_str, "GO Terms": go_terms_desc_str})    
        
            for merged_seq_id in merged_sequences:
                SeqIO.write(merged_sequences[merged_seq_id], out, "fasta")

    df_step3 = pd.DataFrame.from_dict(step3_data).drop_duplicates(subset=["ID", "GO ID"], keep="first")
    df_total = pd.merge(df_total, df_step3, on="ID", how="left")
    df_total["deg.sim"] = df_total["deg.sim"].replace("uncheck", np.nan)
    df_total.loc[(df_total["deg.sim"] > 0) & (df_total["GO Terms"].isnull()), ["GO Terms", "GO ID"]] = "unfit"
    df_total["GO Terms"].fillna("uncheck", inplace=True)
    df_total["GO ID"].fillna("uncheck", inplace=True)
    df_total["deg.sim"].fillna("uncheck", inplace=True)

    ################## STEP4 Perform PFAM annotations on the sequences ############################

    pfam_db_path = os.path.dirname(pfam_db)
    pfam_db_name = os.path.basename(pfam_db)
    if not os.path.exists(f'{pfam_db}.h3f'):
        subprocess.call(f"hmmpress {pfam_db}", shell=True)

    subprocess.call(f"hmmscan --domtblout {step4_file}.domtblout --noali -E 0.01 {pfam_db} {step3_file} > {step4_file}.log", shell=True)
    seq_dict = SeqIO.to_dict(SeqIO.parse(step3_file, "fasta"))

    step4_data = []
    with open(step4_file, "w") as out:
        with open(f"{step4_file}.domtblout", "r") as f:
            domain_dict = {}
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split()
                seq_id = fields[3]
                domain_id = fields[1]
                domain_name = fields[0]
                #print(domain_name)
                start = int(fields[17])
                end = int(fields[18])
                if seq_id not in domain_dict:
                    domain_dict[seq_id] = []
                domain_dict[seq_id].append((domain_id, start, end))
                step4_data.append({"ID": seq_id, "Domain ID": domain_id, "Domain Name": domain_name})

            for seq_id, domains in domain_dict.items():
                # Merge overlapping or close domains
                domains = merge_domains(domains, threshold=30)
                domains.sort(key=lambda x: x[1])
                seq_record = seq_dict[seq_id]
                for i in range(len(domains)):
                    domain_id = domains[i][0]
                    start = domains[i][1]
                    end = domains[i][2]
                    domain_seq = seq_record[start-1:end]
                    out_id = f"{seq_id} {domain_id}"
                    domain_seq.id = out_id
                    domain_seq.description = ""

                seq_record_with_domain = SeqIO.SeqRecord(
                    seq=domain_seq.seq,
                    id=out_id,
                    description=""
                )
                SeqIO.write(seq_record_with_domain, out, "fasta")

    df_step4 = pd.DataFrame.from_dict(step4_data)
    df_step4 = df_step4.groupby('ID').agg({'Domain ID': lambda x: ';'.join(x),
                                        'Domain Name': lambda x: ';'.join(x)}).reset_index()

    df_total = pd.merge(df_total, df_step4, on="ID", how="left")
    df_total.loc[((df_total["GO Terms"] != "unfit") & (df_total["GO Terms"] != "uncheck") & 
                df_total["Domain ID"].isnull() & df_total["Domain Name"].isnull()),
                ["Domain ID", "Domain Name"]] = "unknown"
    df_total["Domain ID"].fillna("uncheck", inplace=True)
    df_total["Domain Name"].fillna("uncheck", inplace=True)

    ################## STEP5 Non-host functional domain identification and retrieval of original pathogenic protein sequences using domain IDs ############################

    # Find non-host domain sequences and save to step5_file_1
    seq_dict = SeqIO.to_dict(SeqIO.parse(step3_file, "fasta"))
    nohost_domain_seqs = []

    # Process the temporary file and save non-host domain sequences to step5_file_1
    step5_data = []
    with open(step5_file_1, "w") as out:
        result = subprocess.run(f"diamond blastp -q {step4_file} -d hostdb --outfmt 6 qseqid sseqid pident evalue --id 90 --subject-cover 80 --query-cover 70 --max-target-seqs 1", shell=True, capture_output=True)
        if result.returncode == 0:
            if result.stdout:
                with open(f"{step5_file_1}.tmp", "w") as f:
                    f.write(result.stdout.decode())
                with open(f"{step5_file_1}.tmp", "r") as f:
                    for line in f:
                        fields = line.strip().split()
                        seq_id = fields[0]
                        subject_id = fields[1]
                        domain_sim = fields[2]
                        evalue = float(fields[3])
                        step5_data.append({"ID": seq_id, "domain.sim": domain_sim})
                        if seq_id in seq_dict and evalue <= 0.01:
                            domain_seq = seq_dict[seq_id]
                            domain_id = seq_id.split("_")[-1]
                            out_id = f"{domain_seq.description} {domain_id}"
                            domain_seq.id = out_id
                            domain_seq.description = ""
                            SeqIO.write(domain_seq, out, "fasta")
                            nohost_domain_seqs.append(domain_seq)
            else:
                with open(step4_file, "r") as f:
                    for record in SeqIO.parse(f, "fasta"):
                        SeqIO.write(record, out, "fasta")
                        step5_data.append({"ID": record.id, "domain.sim": 0})
        else:
            print("Error: diamond blastp failed!")

    # Find corresponding full-length protein sequences and save to step5_file_2
    with open(step5_file_2, "w") as out:
        seq_ids = set()
        with open(step5_file_1, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_ids.add(record.id.split(" ")[0])
        for seq_record in SeqIO.parse(step3_file, "fasta"):
            seq_id = seq_record.id.split(" ")[0]
            if seq_id in seq_ids:
                orig_seq = seq_dict[seq_id]
                out_id = seq_id
                start = orig_seq.seq.find(seq_record.seq)
                end = start + len(seq_record.seq)
                domain_seq = orig_seq[start:end]
                domain_seq.id = out_id
                domain_seq.description = ""
                SeqIO.write(domain_seq, out, "fasta")

    df_step5 = pd.DataFrame.from_dict(step5_data)
    df_total = pd.merge(df_total, df_step5, on="ID", how="left")
    df_total.loc[(df_total["Domain ID"] != "uncheck") & df_total["domain.sim"].isnull(), "domain.sim"] = "0"
    df_total["domain.sim"].fillna("uncheck", inplace=True)

    ################## STEP6 Perform drugability assessment on the sequences ############################

    subprocess.call(f"python spider.py {step5_file_2} output.csv", shell=True)
    df_step6 = pd.read_csv("output.csv", header=None)
    druggable_proteins = df_step6.loc[df_step6.iloc[:, 1] == "Positive", df_step6.columns[0]].tolist()
    with open(step6_file, "w") as f:
        with open(step5_file_2) as seq_file:
            seq_id = None
            seq = ""
            for line in seq_file:
                if line.startswith(">"):
                    if seq_id is not None and seq_id in druggable_proteins:
                        f.write(f">{seq_id}\n{seq}\n")
                    seq_id = line.strip()[1:]
                    seq = ""
                else:
                    seq += line.strip()
            if seq_id is not None and seq_id in druggable_proteins:
                f.write(f">{seq_id}\n{seq}\n")

    df_step6.columns = ["ID", "Druggable", "Predict Score"]
    df_total = pd.merge(df_total, df_step6, on="ID", how="left")
    df_total["Druggable"].fillna("uncheck", inplace=True)
    df_total["Predict Score"].fillna("uncheck", inplace=True)

    # Save the DataFrame to a CSV file
    df_total.to_csv(final_file, index=False)

    # Remove temporary file
    subprocess.call("rm -rf temp", shell=True)
    subprocess.call("rm -rf __pycache__", shell=True)
    subprocess.call("rm output.csv", shell=True)

    print(f"The drug target screening task for the pathogens organism has been completed.")




