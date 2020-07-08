import os
import subprocess
import shutil

from pprint import pprint,pformat

class GATKUtils:

    def __init__(self):
       self.path = "/kb/module/deps/"
       output_dir = "/kb/module/work/tmp/genome/reference/"
       src = "/kb/module/test/genome"

       #self.path = "/home/manish/Desktop/apps/kb_GATK/deps/"
       #output_dir = "/home/manish/Desktop/apps/kb_GATK/test_local/workdir/tmp/genome/reference/"
       #src = "/home/manish/Desktop/apps/kb_GATK/test/genome"

       #assembly_file = os.path.join(src, "reference/NC_008253.fna")
       #fwd_fastq = "/home/manish/Desktop/apps/kb_GATK/test/bt_test_data/reads_1.fq"
       #rev_fastq = "/home/manish/Desktop/apps/kb_GATK/test/bt_test_data/reads_2.fq"
       
       fwd_fastq = "/kb/module/test/bt_test_data/reads_1.fq"
       rev_fastq = "/kb/module/test/bt_test_data/reads_2.fq"
       pass 

    def run_cmd(self, cmd):
        try:
           process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
           stdout, stderr = process.communicate()
           if stdout:
               print ("ret> ", process.returncode)
               print ("OK> output ", stdout)
           if stderr:
               print ("ret> ", process.returncode)
               print ("Error> error ", stderr.strip())

        except OSError as e:
           print ("OSError > ", e.errno)
           print ("OSError > ", e.strerror)
           print ("OSError > ", e.filename)
    
    def build_genome(self, assembly_file ):
        cmd = "bwa index -a bwtsw " + assembly_file
        self.run_cmd(cmd)
       
    def index_assembly(self, assembly_file):
        cmd = "samtools faidx " + assembly_file
        self.run_cmd(cmd)
      
    def generate_sequence_dictionary(self, assembly_file):
       
        output_dict = assembly_file.replace("fa","dict")
        cmd = "java -jar "+ self.path + "picard.jar CreateSequenceDictionary REFERENCE=" + assembly_file +" OUTPUT=" + output_dict
        print(cmd)
        self.run_cmd(cmd)
               
    def mapping_genome(self, ref_genome, rev_fastq, fwd_fastq, output_dir ):
        
        cmd = "bwa mem -t 32 -M -R " + "\"@RG\\tID:sample_1\\tLB:sample_1\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:sample_1\" " + ref_genome +" "+ rev_fastq + " " + fwd_fastq + " > " + output_dir + "aligned_reads.sam"
        self.run_cmd(cmd)
    
    def duplicate_marking(self, output_dir, sam_file):
        cmd = "java -jar "+ self.path + "picard.jar SortSam  INPUT= " + sam_file + "  OUTPUT=" + output_dir + "aligned_reads.bam  SORT_ORDER=coordinate"
        #cmd = "java -jar "+ self.path + "picard.jar SortSam  INPUT= " + output_dir + "aligned_reads.sam   OUTPUT=" + output_dir + "aligned_reads.bam  SORT_ORDER=coordinate"
        self.run_cmd(cmd)
       
    def sort_bam_index(self, output_dir):
        cmd = "samtools index " + output_dir + "aligned_reads.bam"
        self.run_cmd(cmd)
        
    def collect_alignment_and_insert_size_metrics(self, assembly_file, output_dir):
        cmd1 = "java -jar "+ self.path + "picard.jar CollectAlignmentSummaryMetrics R="+ assembly_file  +" I=" + output_dir + "aligned_reads.bam O=" + output_dir + "alignment_metrics.txt"
        self.run_cmd(cmd1)
        cmd2 = "java -jar "+ self.path + "picard.jar CollectInsertSizeMetrics INPUT=" + output_dir + "aligned_reads.bam OUTPUT=" + output_dir + "insert_metrics.txt HISTOGRAM_FILE=" + output_dir + "insert_size_histogram.pdf"
        self.run_cmd(cmd2)
        cmd3 = "samtools depth -a " + output_dir + "aligned_reads.bam > " + output_dir + "depth_out.txt"
        self.run_cmd(cmd3)

    def variant_calling(self, assembly_file, output_dir):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar HaplotypeCaller -R "+ assembly_file  + " -I " + output_dir + "aligned_reads.bam -O " + output_dir + "raw_variants.vcf"
        self.run_cmd(cmd)

    def extract_variants(self, assembly_file, output_dir ):
        cmd1 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants -R " + assembly_file  + " -V " + output_dir + "raw_variants.vcf --select-type SNP -O " + output_dir + "raw_snps.vcf"
        self.run_cmd(cmd1)
        cmd2 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants -R " + assembly_file  + " -V " + output_dir + "raw_variants.vcf --select-type INDEL -O " + output_dir + "raw_indels.vcf"
        self.run_cmd(cmd2)

    def filter_SNPs(self, assembly_file, output_file, output_dir, params):
        print(params)
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration -R " + assembly_file  + " -V " + output_dir + "raw_snps.vcf -O "  + output_dir +  output_file +" -filter-name 'QD_filter' -filter 'QD < " + params['snp_qd_filter'] + "' -filter-name 'FS_filter' -filter 'FS > " + params['snp_fs_filter'] + "' -filter-name 'MQ_filter' -filter 'MQ < " + params['snp_mq_filter'] + "' -filter-name 'SOR_filter' -filter 'SOR > " + params['snp_sor_filter'] + "' -filter-name 'MQRankSum_filter' -filter 'MQRankSum < " + params['snp_mqrankSum_filter'] + "' -filter-name 'ReadPosRankSum_filter' -filter 'ReadPosRankSum < " + params['snp_readposranksum_filter'] + "'"
        print(cmd)
        self.run_cmd(cmd)

    def filter_Indels(self, assembly_file, output_file, output_dir, params):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration -R " + assembly_file  + " -V " + output_dir + "raw_indels.vcf -O " + output_dir + output_file +" -filter-name 'QD_filter' -filter 'QD < " + params['indel_qd_filter'] + "' -filter-name 'FS_filter' -filter 'FS > " + params['indel_fs_filter'] + "' -filter-name 'SOR_filter' -filter 'SOR > " + params['indel_sor_filter'] + "'"
        self.run_cmd(cmd)

    def exclude_filtered_variants(self, output_dir):
        cmd1 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants --exclude-filtered -V "  + output_dir +"filtered_snps.vcf -O "  + output_dir +"bqsr_snps.vcf"
        self.run_cmd(cmd1)

        cmd2 = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar SelectVariants --exclude-filtered -V "  + output_dir + "filtered_indels.vcf -O "  + output_dir +"bqsr_indels.vcf"
        self.run_cmd(cmd2)

    def base_quality_score_recalibration(self, assembly_file, data_table, output_dir):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar BaseRecalibrator -R " + assembly_file  + " -I "  + output_dir + "aligned_reads.bam --known-sites "  + output_dir + "bqsr_snps.vcf --known-sites "  + output_dir + "bqsr_indels.vcf -O "  + output_dir + data_table
        self.run_cmd(cmd)

    def apply_BQSR(self, assembly_file, data_table, output_dir):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar ApplyBQSR -R " + assembly_file  + "  -I "  + output_dir +"aligned_reads.bam -bqsr "+ output_dir + data_table + " -O "  + output_dir +"recal_reads.bam"
        self.run_cmd(cmd)

    def analyze_covariates(self, output_dir):
        cmd = "java -jar "+ self.path + "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar AnalyzeCovariates -before "  + output_dir +"recal_data.table -after "  + output_dir +"post_recal_data.table -plots "  + output_dir +"recalibration_plots.pdf"
        self.run_cmd(cmd)
'''
if __name__ == "__main__":
   
  gu = GATKUtils() 
  src = "/home/manish/Desktop/apps/kb_GATK/test/genome"
  output_dir = "/home/manish/Desktop/apps/kb_GATK/test_local/workdir/tmp"
  #output_dir = os.path.join(self.shared_folder, str(uuid.uuid4()))
  #os.mkdir(output_dir)
  #dest = shutil.copytree(src, os.path.join(output_dir, "genome"))

  assembly_file = os.path.join(src, "reference/test.fna")
       
  #gu.build_genome(assembly_file)

  #gu.index_assembly(assembly_file)

  #gu.generate_sequence_dictionary(assembly_file)

  fwd_fastq = "/home/manish/Desktop/apps/kb_GATK/test/bt_test_data/reads_1.fq"
  rev_fastq = "/home/manish/Desktop/apps/kb_GATK/test/bt_test_data/reads_2.fq"

  #gu.mapping_genome(assembly_file, fwd_fastq, rev_fastq )
 
 
  #gu.duplicate_marking()
   
  #gu.sort_bam_index()
  
  #gu.collect_alignment_and_insert_size_metrics(assembly_file)

  gu.analyze_covariates()
  
  #gu.variant_calling(assembly_file)
 
  #gu.extract_variants(assembly_file)
   
  #gu.filter_SNPs(assembly_file, "filtered_snps.vcf")
   
  #gu.filter_Indels(assembly_file, "filtered_indels.vcf")
    
  #gu.exclude_filtered_variants()
   
  #gu.base_quality_score_recalibration(assembly_file, "recal_data.table")
  
  #gu.apply_BQSR(assembly_file, "recal_data.table")
  
  #gu.base_quality_score_recalibration(assembly_file, "post_recal_data.table")
    
  #gu.apply_BQSR(assembly_file,  "post_recal_data.table")
   
  #gu.filter_SNPs(assembly_file, "filtered_snps_final.vcf")
   
  #gu.filter_Indels(assembly_file, "filtered_indels_final.vcf")
  
'''

