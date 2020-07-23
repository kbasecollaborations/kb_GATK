import subprocess


class GATKUtils:

    def __init__(self):
        self.path = "/kb/module/deps/"
        pass

    def run_cmd(self, cmd):
        try:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if stdout:
                print("ret> ", process.returncode)
                print("OK> output ", stdout)
            if stderr:
                print("ret> ", process.returncode)
                print("Error> error ", stderr.strip())

        except OSError as e:
            print("OSError > ", e.errno)
            print("OSError > ", e.strerror)
            print("OSError > ", e.filename)

    def build_genome(self, assembly_file):
        cmd = "bwa index -a bwtsw " + assembly_file
        self.run_cmd(cmd)

    def index_assembly(self, assembly_file):
        cmd = "samtools faidx " + assembly_file
        self.run_cmd(cmd)

    def generate_sequence_dictionary(self, assembly_file):
        output_dict = assembly_file.replace("fa", "dict")
        cmd = "java -jar " + self.path + "picard.jar CreateSequenceDictionary REFERENCE=" + assembly_file + " OUTPUT=" + output_dict
        self.run_cmd(cmd)

    def mapping_genome(self, ref_genome, rev_fastq, fwd_fastq, output_dir):
        cmd = "bwa mem -t 32 -M -R " + "\"@RG\\tID:sample_1\\tLB:sample_1\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:sample_1\" " + ref_genome + " " + rev_fastq + " " + fwd_fastq + " > " + output_dir + "aligned_reads.sam"
        self.run_cmd(cmd)

    def duplicate_marking(self, output_dir, sam_file):
        cmd = "java -jar " + self.path + "picard.jar SortSam  INPUT= " + sam_file + "  OUTPUT=" + output_dir + "aligned_reads.bam  SORT_ORDER=coordinate"
        self.run_cmd(cmd)

    def sort_bam_index(self, output_dir):
        cmd = "samtools index " + output_dir + "aligned_reads.bam"
        self.run_cmd(cmd)

    def collect_alignment_and_insert_size_metrics(self, assembly_file, output_dir):
        cmd1 = "java -jar " + self.path + "picard.jar CollectAlignmentSummaryMetrics R=" + assembly_file + " I=" + output_dir + "aligned_reads.bam O=" + output_dir + "alignment_metrics.txt"
        self.run_cmd(cmd1)
        cmd2 = "java -jar " + self.path + "picard.jar CollectInsertSizeMetrics INPUT=" + output_dir + "aligned_reads.bam OUTPUT=" + output_dir + "insert_metrics.txt HISTOGRAM_FILE=" + output_dir + "insert_size_histogram.pdf"
        self.run_cmd(cmd2)
        cmd3 = "samtools depth -a " + output_dir + "aligned_reads.bam > " + output_dir + "depth_out.txt"
        self.run_cmd(cmd3)

    def variant_calling(self, assembly_file, output_dir):
        cmd = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar HaplotypeCaller -R " + assembly_file + " -I " + output_dir + "aligned_reads.bam -O " + output_dir + "raw_variants.vcf"
        self.run_cmd(cmd)

    def extract_variants(self, assembly_file, output_dir):
        cmd1 = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants -R " + assembly_file + " -V " + output_dir + "raw_variants.vcf --select-type SNP -O " + output_dir + "raw_snps.vcf"
        self.run_cmd(cmd1)
        cmd2 = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants -R " + assembly_file + " -V " + output_dir + "raw_variants.vcf --select-type INDEL -O " + output_dir + "raw_indels.vcf"
        self.run_cmd(cmd2)

    def filter_SNPs(self, assembly_file, output_file, output_dir, params):
        # Todo: need to remove hardcoded filters

        '''
        params['snp_filter']['snp_qd_filter'] = '2.0'
        params['snp_filter']['snp_fs_filter'] = '60.0'
        params['snp_filter']['snp_mq_filter'] = '40.0'
        params['snp_filter']['snp_sor_filter'] = '4.0'
        params['snp_filter']['snp_mqrankSum_filter'] = '-12.5'
        params['snp_filter']['snp_readposranksum_filter'] = '-8.0'
        '''

        cmd = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration -R " + assembly_file + " -V " + output_dir + "raw_snps.vcf -O " + output_dir + output_file + " -filter-name 'QD_filter' -filter 'QD < " + params['snp_filter']['snp_qd_filter'] + "' -filter-name 'FS_filter' -filter 'FS > " + params['snp_filter']['snp_fs_filter'] + "' -filter-name 'MQ_filter' -filter 'MQ < " + params['snp_filter']['snp_mq_filter'] + "' -filter-name 'SOR_filter' -filter 'SOR > " + params['snp_filter']['snp_sor_filter'] + "' -filter-name 'MQRankSum_filter' -filter 'MQRankSum < " + params['snp_filter']['snp_mqrankSum_filter'] + "' -filter-name 'ReadPosRankSum_filter' -filter 'ReadPosRankSum < " + params['snp_filter']['snp_readposranksum_filter'] + "'"
        self.run_cmd(cmd)

    def filter_Indels(self, assembly_file, output_file, output_dir, params):
        # Todo: need to removed hardcoded filters

        '''
        params['indel_filter']['indel_qd_filter'] = '2.0'
        params['indel_filter']['indel_fs_filter'] = '200.0'
        params['indel_filter']['indel_mq_filter'] = '40.0'
        params['indel_filter']['indel_sor_filter'] = '10.0'
        params['indel_filter']['indel_mqrankSum_filter'] = '-12.5'
        params['indel_filter']['indel_readposranksum_filter'] = '-8.0'
        '''

        cmd = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration -R " + assembly_file + " -V " + output_dir + "raw_indels.vcf -O " + output_dir + output_file + " -filter-name 'QD_filter' -filter 'QD < " + params['indel_filter']['indel_qd_filter'] + "' -filter-name 'FS_filter' -filter 'FS > " + params['indel_filter']['indel_fs_filter'] + "' -filter-name 'SOR_filter' -filter 'SOR > " + params['indel_filter']['indel_sor_filter'] + "'"
        self.run_cmd(cmd)

    def exclude_filtered_variants(self, output_dir):
        cmd1 = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants --exclude-filtered -V " + output_dir + "filtered_snps.vcf -O " + output_dir + "bqsr_snps.vcf"
        self.run_cmd(cmd1)

        cmd2 = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants --exclude-filtered -V " + output_dir + "filtered_indels.vcf -O " + output_dir + "bqsr_indels.vcf"
        self.run_cmd(cmd2)

    def base_quality_score_recalibration(self, assembly_file, data_table, output_dir):
        cmd = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar BaseRecalibrator -R " + assembly_file + " -I " + output_dir + "aligned_reads.bam --known-sites " + output_dir + "bqsr_snps.vcf --known-sites " + output_dir + "bqsr_indels.vcf -O " + output_dir + data_table
        self.run_cmd(cmd)

    def apply_BQSR(self, assembly_file, data_table, output_dir):
        cmd = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar ApplyBQSR -R " + assembly_file + "  -I " + output_dir + "aligned_reads.bam -bqsr " + output_dir + data_table + " -O " + output_dir + "recal_reads.bam"
        self.run_cmd(cmd)

    def analyze_covariates(self, output_dir):
        cmd = "java -jar " + self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar AnalyzeCovariates -before " + output_dir + "recal_data.table -after " + output_dir + "post_recal_data.table -plots " + output_dir + "recalibration_plots.pdf"
        self.run_cmd(cmd)
