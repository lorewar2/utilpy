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
VCF_LOC = "/data1/phasstphase_test/potato/potato_deep_phased.vcf.gz"
INTERMEDIATE_LOC = "./intermediate/"
HERA_UNIQUE_LOC = "./intermediate/hera_unique.fa.ktab"
STIEG_UNIQUE_LOC = "./intermediate/stieg_unique.fa.ktab"
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
    #util.find_which_parent_contain_kstring(0, k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, TABEX_LOC, INTERMEDIATE_LOC, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC)
    #util.thread_runner_kmer_search(THREAD_NUMBER, k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, TABEX_LOC, INTERMEDIATE_LOC, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC)
    #util.look_for_stieg_ref(0, k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, TABEX_LOC, INTERMEDIATE_LOC, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC)
    #util.find_specific_phaseblock_kmer(K, VCF_LOC, HERA1_REF_LOC, 156, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC) #41039600 #156 #3224262
    #util.make_fasta_file_for_each_phase_block_haplotype(K, VCF_LOC, HERA1_REF_LOC, "/data1/phasstphase_test/potato/kmc_run/phaseblockfasta")
    #util.make_kmc_files_and_dump("/data1/phasstphase_test/potato/kmc_run/phaseblockfasta")
    #util.make_result_file_from_dump ("/data1/phasstphase_test/potato/kmc_run/phaseblockfasta")
    util.save_alt_and_ref_kmers_in_files(21, "/data1/GiaB_benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer_fix.resorted.vcf", "/data1/GiaB_benchmark/GRCh38.fa")
    return

if __name__ == "__main__":
    main()