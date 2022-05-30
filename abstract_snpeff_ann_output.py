import re, os
from shutil import copyfile

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
fout = open("k-per-gene_variant_anns.tsv",'w')
fout.write(header)
suffix = '.snpeff.vcf'
for file in os.listdir():
#file = "COVC14272_S10_L001_variants.ann.vcf"
    if file.endswith(suffix):
        # print(file)
        cov_name = file.split(".")[0].split("_")[0]
        sam_name = cov_name#"_".join(cov_name)
        with open(file) as f:
            genes = {'orf1ab':[], 'ORF1a':[], 'S':[], 'ORF3a':[], 'ORF3b':[]
            , 'E':[], 'M':[], 'ORF6':[], 'ORF7a':[], 'ORF7b':[], 'ORF8':[]
            , 'N':[], 'ORF9a':[], 'ORF9b':[], 'ORF10':[]}#, 'ORF6&ORF7a': []}
            for line in f:
                if not re.match("#",line):
                    line = re.split("\t", line)
                    ann = line[7]
                    gene = ann.split("|")[3]
                    pos = "g." + (str(line[1]))
                    HGVS_c = ann.split("|")[9]
                    HGVS_p = replace(ann.split("|")[10], aa_code)
                    prot_mut = [HGVS_p][0].lstrip('p.')
                    if prot_mut == '': prot_mut = "NC"
                    try:
                        genes[gene].append(prot_mut)
                    except:
                        if KeyError:
                            print(f'\nWARNING: A unique gene feature ({gene}) was noted in sample {sam_name}\n')
                    pass
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
            # print(out_line)
            ncol = 15
            if out_line != sam_name + '\t' + str(var_count) + ncol*'\t' + '\n':
                fout.write(out_line)
            else: 
                print(f"Sample {sam_name} didn't have any variants called")
                fout.write(out_line)
    elif "k-per-gene_variant_anns.tsv" not in file:
        print(f"File '{file}' is not a '{suffix}' file type: It was skipped...")
fout.close()
copyfile("./k-per-gene_variant_anns.tsv", "../var/k-per-gene_variant_anns.tsv")