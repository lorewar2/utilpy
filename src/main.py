import util

FASTK_LOC = "/home/mweerakoon/fastk/FASTK/FastK"
TABEX_LOC = "/home/mweerakoon/fastk/FASTK/Tabex"
HERA1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera1/solTubHeraHap1.fa"
HERA2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera2/solTubHeraHap2.fa"
STIEG1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg1/solTubStiegHap1.fa"
STIEG2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg2/solTubStiegHap2.fa"
VCF_LOC = "/data1/phasstphase_test/potato/hifi/potato_deep.vcf"
INTERMEDIATE_LOC = "./intermediate/"
K = 30

def main():
    # make intermediate files for hera12 stieg12
    #util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, HERA1_REF_LOC)
    #util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, HERA2_REF_LOC)
    #util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, STIEG1_REF_LOC)
    #util.run_fastk_make_intermediate_files(K, FASTK_LOC, INTERMEDIATE_LOC, STIEG2_REF_LOC)
    # test
    # search the vcf and find a kmer

    # look for the kmer in all 4 parents
    #util.search_for_kmer_in_intermediate(TABEX_LOC, INTERMEDIATE_LOC, HERA1_REF_LOC, "aaaaaaatatttagtggtgataaattttct")
    k_string_vec = util.open_vcf_and_get_k_mer(K, VCF_LOC, HERA1_REF_LOC)
    util.find_which_parent_contain_kstring(k_string_vec, TABEX_LOC, INTERMEDIATE_LOC, [HERA1_REF_LOC, HERA2_REF_LOC, STIEG1_REF_LOC, STIEG2_REF_LOC])
    return

if __name__ == "__main__":
    main()