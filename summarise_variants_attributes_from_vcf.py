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

def get_var_attr(x):
    prot_mut = '';alt_freq = '';alt_qual = '';g_coverage = ''; g_pos = ''
    prot_mut += x[0]
    alt_freq += str(x[1])
    alt_qual += str(x[2])
    g_coverage += str(x[3])
    g_pos += str(x[4])
    return prot_mut, alt_freq, alt_qual, g_coverage, g_pos

header = 'sample_name\tnum_vars\tORF1ab\tORF1a\tS\tORF3a\tORF3b\tE\tM\tORF6\tORF7a\tORF7b\tORF8\tN\tORF9a\tORF9b\tORF10\n'
fout = open("per-gene_variants_QC-attr.tsv",'w')
fout.write(header)
suffix = '.snpeff.vcf'
for file in os.listdir():
#file = "COVC14272_S10_L001_variants.ann.vcf"
    if file.endswith(suffix):
        # print(file)
        cov_name = file.split(".")[0]
        sam_name = cov_name#"_".join(cov_name)
        with open(file) as f:
            genes = {'orf1ab':[], 'ORF1a':[], 'S':[], 'ORF3a':[], 'ORF3b':[]
            , 'E':[], 'M':[], 'ORF6':[], 'ORF7a':[], 'ORF7b':[], 'ORF8':[]
            , 'N':[], 'ORF9a':[], 'ORF9b':[], 'ORF10':[]}
            for line in f:
                if not re.match("#",line) and line != '':
                    line = re.split("\t", line)
                    ann = line[7]
                    if line[8].split(':') == ['GT', 'GQ', 'PS', 'UG', 'UQ']:
                        coverage = line[7].split(';')[0].split('=')[1]
                        alt_freq = round((int(line[7].split(';')[1].split(',')[1]) / int(coverage))*100)
                        alt_qual = round(float(line[7].split(';')[6].split('=')[1]))
                    else:
                        coverage = int(line[9].split(':')[1]) + int(line[9].split(':')[4])
                        alt_freq = round(float(line[9].split(':')[-1]) * 100)
                        alt_qual = line[9].split(':')[6]
                    gene = ann.split("|")[3]
                    pos = "g." + (str(line[1]))
                    HGVS_c = ann.split("|")[9]
                    HGVS_p = replace(ann.split("|")[10], aa_code)
                    prot_mut = [HGVS_p][0].lstrip('p.')
                    if prot_mut == '': prot_mut = "NC"
                    vcf_tup = (prot_mut, alt_freq, alt_qual, coverage, pos)
                    genes[gene].append(vcf_tup)
            variants = ''; altfreq = ''; altqual = ''; gencov = ''; genpos = ''
            var_count = 0
            for gene in genes:
                var_attrs = map(get_var_attr, genes[gene])
                vart = ''; altf = ''; altq = ''; genc = ''; genp = ''
                for var_attr in var_attrs:
                    var_count += 1
                    pm,af,aq,gc,gp = var_attr
                    vart += f', {pm}'
                    altf += f', {af}'#frequency
                    altq += f', {aq}'#quality
                    genc += f', {gc}'#coverage
                    genp += f', {gp}'#position
                variants = variants + '\t' + vart.lstrip(', ')
                altfreq = altfreq + '\t' + altf.lstrip(', ')
                altqual = altqual + '\t' + altq.lstrip(', ')
                gencov = gencov + '\t' + genc.lstrip(', ')
                genpos = genpos + '\t' + genp.lstrip(', ')

            out_line1 = sam_name + '\t' + str(var_count) + variants + '\n'
            out_line2 = f'af_{sam_name}' + '\t' + str(var_count) + altfreq + '\n'
            out_line3 = f'aq_{sam_name}' + '\t' + str(var_count) + altqual + '\n'
            out_line4 = f'gc_{sam_name}' + '\t' + str(var_count) + gencov + '\n'
            out_line5 = f'gp_{sam_name}' + '\t' + str(var_count) + genpos + '\n'
            # print(out_line)
            if out_line1 != sam_name + '\t' + str(var_count) + 15*'\t' + '\n':
                fout.write(out_line1)
                fout.write(out_line2)
                fout.write(out_line3)
                fout.write(out_line4)
                fout.write(out_line5)
            else: 
                print(f"Sample {sam_name} didn't have any variants called")
                fout.write(out_line1)
    else: print(f"File '{file}' is not a VCF file type: It was skipped...")
fout.close()