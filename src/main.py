import util

FASTK_LOC = "/home/mweerakoon/fastk/FASTK/FastK"
TABEX_LOC = "/home/mweerakoon/fastk/FASTK/Tabex"
LOGEX_LOC = "/home/mweerakoon/fastk/FASTK/Logex"
HERA1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera1/solTubHeraHap1.fa"
HERA2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera2/solTubHeraHap2.fa"
STIEG1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg1/solTubStiegHap1.fa"
STIEG2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg2/solTubStiegHap2.fa"
HERA_READ_LOC = "/data1/phasstphase_test/potato/hera/hera_merged.fastq"
STIEG_READ_LOC = "/data1/phasstphase_test/potato/stieglitz/stieg_merged.fastq"
DEEPVARIANT_LOC = "/data1/phasstphase_test/potato/hifi/potato_deep.g.vcf.gz"
VCF_LOC = "/data1/phasstphase_test/potato/potato_deep_phased.vcf.gz"
INTERMEDIATE_LOC = "./intermediate/"
HERA_UNIQUE_LOC = "./intermediate/hera_unique.fa.ktab"
STIEG_UNIQUE_LOC = "./intermediate/stieg_unique.fa.ktab"
GROUND1_SAVE_LOC = "intermediate/hap1_ground.fa"
GROUND2_SAVE_LOC = "intermediate/hap2_ground.fa"
THREAD_NUMBER = 1
K = 21

def main():
    # make intermediate files for merged hera and merged stieg
    #util.run_fastk_make_intermediate_files_from_reads(K, FASTK_LOC, INTERMEDIATE_LOC, HERA_READ_LOC, 6)
    #util.run_fastk_make_intermediate_files_from_reads(K, FASTK_LOC, INTERMEDIATE_LOC, STIEG_READ_LOC, 5)
    # make unique kmer files
    #util.unique_kmers_for_parent_from_intermediates_from_reads(LOGEX_LOC, INTERMEDIATE_LOC, HERA_READ_LOC, STIEG_READ_LOC)
    # search the vcf and find a kmer
    #k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks = util.open_vcf_and_get_k_mer(K, VCF_LOC, HERA1_REF_LOC)
    # look for the kmer in all 4 parents
    util.find_specific_phaseblock_kmer(K, VCF_LOC, HERA1_REF_LOC, 3224262, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC)
    #util.find_which_parent_contain_kstring(0, k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, TABEX_LOC, INTERMEDIATE_LOC, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC)
    return

if __name__ == "__main__":
    main()