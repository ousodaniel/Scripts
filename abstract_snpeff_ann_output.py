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

header = 'sample_name\tannotation\n'
fout = open("All_variant_anns.tsv",'w')
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
                    cp = tuple([pos, HGVS_c, HGVS_p])
                    genes[gene].append(cp)
            ann_str = ""
            for gene in genes:
                if len(genes[gene]) > 0:
                    x = ('{}: {}'.format(gene, str(genes[gene])
                                 .replace("[", "")
                                 .replace("]", "")
                                 .replace("'", ""))
                        )
                    ann_str = ann_str + x + '\n\t'
                else: pass
            out_line = sam_name + '\t' + ann_str + '\n'
            if out_line == sam_name + '\t\n':
                pass
            else: fout.write(out_line)
fout.close()
