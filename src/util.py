import os
import vcf
import pyfaidx

def run_fastk_make_intermediate_files(k, fast_k_loc, intermediate_loc, ref_loc):
    # get the file name from path
    file_name = ref_loc.split("/")[-1]
    intermediate_path = "{}{}".format(intermediate_loc, file_name)
    command_to_run = "{} -k{} -t1 -N\"{}\" -T64 {}".format(fast_k_loc, k, intermediate_path, ref_loc)
    # run the command
    print(command_to_run)
    print(os.popen(command_to_run).read())
    return

def open_vcf_and_get_k_mer(k, vcf_loc, ref_loc):
    ref_alt_kmer_list = []
    variant_reader = vcf.Reader(filename = vcf_loc)
    ref_fasta = pyfaidx.Fasta(ref_loc)
    if k % 2 == 0:
        k_first_half_length = k // 2
        k_second_half_length = k // 2
    else:
        k_first_half_length = (k // 2) + 1
        k_second_half_length = k // 2
    for record in variant_reader:
        if len(record.alleles) > 2:
            continue
        ref = record.alleles[0]
        alt = record.alleles[1]
        kmer_first_half = ref_fasta[record.CHROM][record.POS - k_first_half_length : record.POS]
        kmer_second_half_ref = ref_fasta[record.CHROM][record.POS + len(ref) : record.POS + k_second_half_length]
        kmer_second_half_alt = ref_fasta[record.CHROM][record.POS + len(ref) : record.POS + len(ref) - len(alt) + k_second_half_length]
        if ((len(kmer_first_half) + len(kmer_second_half_ref) + len(ref)) == k) and ((len(kmer_first_half) + len(kmer_second_half_alt) + len(alt)) == k):
            ref_kmer = (kmer_first_half + ref + kmer_second_half_ref).lower()
            alt_kmer = (kmer_first_half + alt + kmer_second_half_alt).lower()
            ref_alt_kmer_list.append((ref_kmer, alt_kmer))
            print(ref_kmer + "appended")
    return ref_alt_kmer_list

def search_for_kmer_in_intermediate(tabex_loc, intermediate_loc, ref_loc, k_string):
    file_name = ref_loc.split("/")[-1]
    ktab_path = "{}{}.ktab".format(intermediate_loc, file_name)
    command_to_run = "{} {} {} {}".format(tabex_loc, ktab_path, ref_loc, k_string)
    print(command_to_run)
    output = os.popen(command_to_run).read()
    if output.find("Not Found") == -1:
        print("KMER DOESNT EXIST")
    else:
        print("KMER EXIST")
    return

#./FastK -k30 -N"./../" -T64 /data1/phasstphase_test/potato/reference/solTubHeraHap1.fa
# aaaaaaatatttagtggtgataaattttct