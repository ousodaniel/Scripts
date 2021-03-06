import re, os

aa_code = {
    "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F",
    "Gly": "G", "His": "H", "Ile": "I", "Lys": "K", "Leu": "L",
    "Met": "M", "Asn": "N", "Pro": "P", "Gln": "Q", "Arg": "R",
    "Ser": "S", "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y",
    "ter": "*"
    }

def replace(string, substitutions):

    substrings = sorted(substitutions, key=len, reverse=True)
    regex = re.compile('|'.join(map(re.escape, substrings)))
    return regex.sub(lambda match: substitutions[match.group(0)], string)

header = 'sample_name\tnum_vars\tORF1ab\tORF1a\tS\tORF3a\tORF3b\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF9a\tORF9b\tORF10\n'
fout = open("per-gene_variant_anns.tsv",'w')
fout.write(header)
for file in os.listdir():
#file = "COVC14272_S10_L001_variants.ann.vcf"
    if file.endswith('.ann.vcf'):
        cov_name = file.split("_")[0:2]
        sam_name = "_".join(cov_name)
        with open(file) as f:
            genes = {'ORF1ab':[], 'ORF1a':[], 'S':[], 'ORF3a':[], 'ORF3b':[]
            , 'E':[], 'M':[], 'ORF6':[], 'ORF7a':[], 'ORF7b':[], 'ORF8':[]
            , 'N':[], 'ORF9a':[], 'ORF9b':[], 'ORF10':[]}
            for line in f:
                if not re.match("#",line):
                    line = re.split("\t", line)
                    ann = line[7]
                    gene = ann.split("|")[3]
                    pos = "g." + (str(line[1]))
                    HGVS_c = ann.split("|")[9]
                    HGVS_p = replace(ann.split("|")[10], aa_code)
                    prot_mut = [HGVS_p][0].lstrip('p.')
                    genes[gene].append(prot_mut)
            ann_str = ""
            var_count = 0
            for gene in genes:
                for var in genes[gene]:
                    if len(var) != 0:
                        var_count += 1
                    else: pass
                _ = ('{}'.format(str(genes[gene])
                             .replace("[", "")
                             .replace("]", "")
                             .replace("'", "")
                             .lstrip(', '))
                    )
                ann_str = ann_str + '\t' + _
            out_line = sam_name + '\t' + str(var_count) + ann_str + '\n'
            if out_line != sam_name + '\t' + str(var_count) + 15*'\t' + '\n':
                fout.write(out_line)
            else: print(f"Sample {sam_name} didn't have any variants callled")
fout.close()