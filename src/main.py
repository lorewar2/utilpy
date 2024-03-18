import util

FASTK_LOC = "/home/mweerakoon/fastk/FASTK/FastK"
TABEX_LOC = "/home/mweerakoon/fastk/FASTK/Tabex"
LOGEX_LOC = "/home/mweerakoon/fastk/FASTK/Logex"
HERA1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera1/solTubHeraHap1.fa"
HERA2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera2/solTubHeraHap2.fa"
STIEG1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg1/solTubStiegHap1.fa"
STIEG2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg2/solTubStiegHap2.fa"
VCF_LOC = "/data1/phasstphase_test/potato/potato_deep_phased.vcf.gz"
INTERMEDIATE_LOC = "./intermediate/"
K = 40

def main():
    # make intermediate files for hera12 stieg12
    util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, HERA1_REF_LOC)
    util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, HERA2_REF_LOC)
    util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, STIEG1_REF_LOC)
    util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, STIEG2_REF_LOC)
    # make unique kmer files
    util.unique_kmers_for_parent_from_intermediates(LOGEX_LOC, INTERMEDIATE_LOC, [HERA1_REF_LOC, HERA2_REF_LOC, STIEG1_REF_LOC, STIEG2_REF_LOC])
    # search the vcf and find a kmer
    #k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks = util.open_vcf_and_get_k_mer(K, VCF_LOC, HERA1_REF_LOC)
    # look for the kmer in all 4 parents
    #util.find_which_parent_contain_kstring(k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, TABEX_LOC, INTERMEDIATE_LOC, [HERA1_REF_LOC, HERA2_REF_LOC, STIEG1_REF_LOC, STIEG2_REF_LOC])
    return

if __name__ == "__main__":
    main()