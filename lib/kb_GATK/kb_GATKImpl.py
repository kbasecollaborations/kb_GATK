# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
from installed_clients.VariationUtilClient import VariationUtil
from installed_clients.KBaseReportClient import KBaseReport
from kb_GATK.Utils.GATKUtils import GATKUtils
from kb_GATK.Utils.DownloadAlignmentUtils import DownloadAlignmentUtils
from installed_clients.WorkspaceClient import Workspace

#END_HEADER


class kb_GATK:
    '''
    Module Name:
    kb_GATK

    Module Description:
    A KBase module: kb_GATK
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbasecollaborations/kb_GATK.git"
    GIT_COMMIT_HASH = "5e6e4bdca9a7749bba0abab081736c56007212ed"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.ws_url = config['workspace-url']
        self.wsc = Workspace(self.ws_url)
        self.gu = GATKUtils()
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.vu = VariationUtil(self.callback_url)
        self.du = DownloadAlignmentUtils(self.callback_url)
        #END_CONSTRUCTOR
        pass

    def run_kb_GATK(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_GATK
        print(params)
        source_ref = params['alignment_ref']
        alignment_out = self.du.downloadreadalignment(source_ref, params, self.callback_url)
        sam_file = os.path.join(alignment_out['destination_dir'], "reads_alignment.sam")

        '''
        #Todo Reading sample set and sample strains information
        '''

        '''
        command.extend(["-filter-name", "\"QD_filter\"", "-filter", "\"QD", "<", params['snp_filter']['snp_qd_filter'] + "\""])
        command.extend(["-filter-name", "\"FS_filter\"", "-filter", "\"FS", "<", params['snp_filter']['snp_fs_filter'] + "\""])
        command.extend(["-filter-name", "\"MQ_filter\"", "-filter", "\"MQ", "<", params['snp_filter']['snp_mq_filter'] + "\""])
        command.extend(["-filter-name", "\"SOR_filter\"", "-filter", "\"SOR", "<", params['snp_filter']['snp_sor_filter'] + "\""])
        command.extend(["-filter-name", "\"MQRankSum_filter\"", "-filter", "\"MQRankSum", "<", params['snp_filter']['snp_mqrankSum_filter'] + "\""])
        command.extend(["-filter-name", "\"ReadPosRankSum_filter\"", "-filter", "\"ReadPosRankSum", "<", params['snp_filter']['snp_readposranksum_filter'] + "\""])
        '''

        output_dir = os.path.join(self.shared_folder, str(uuid.uuid4()))
        os.mkdir(output_dir)

        genome_or_assembly_ref = params['assembly_or_genome_ref']
        obj_type = self.wsc.get_object_info3({
            'objects':[{
                'ref': genome_or_assembly_ref
                      }]})['infos'][0][2]
        if ('KBaseGenomes.Genome' in obj_type):
            genome_ref = genome_or_assembly_ref
            subset = self.wsc.get_object_subset([{
                    'included': ['/assembly_ref'],
                    'ref': genome_ref
                }])
            assembly_ref = subset[0]['data']['assembly_ref']
        elif ('KBaseGenomeAnnotations.Assembly' in obj_type):
            assembly_ref = genome_or_assembly_ref
        else:
            raise ValueError(obj_type + ' is not the right input for this method. '
                                      + 'Valid input include KBaseGenomes.Genome or '
                                      + 'KBaseGenomeAnnotations.Assembly ')       

        assembly_file = self.du.download_genome(assembly_ref, output_dir)['path']

        #output_dir = output_dir + "/"

        #Todo: check time for building index file or donwload from cache. 
        #Todo: To discuss about cache_id to be used.
        #Todo: In case of copying genome, find the way of finding original genome (ref id) for getting original cache id.

        self.gu.build_genome(assembly_file)
        self.gu.index_assembly(assembly_file)
        self.gu.generate_sequence_dictionary(assembly_file)
        self.gu.duplicate_marking(output_dir, sam_file)
        self.gu.sort_bam_index(output_dir)
        self.gu.collect_alignment_and_insert_size_metrics(assembly_file, output_dir)
        #self.gu.analyze_covariates(output_dir)

        #Todo: avoid writing intermediate fies to save space and time I/O. 
        self.gu.variant_calling(assembly_file, output_dir)
        self.gu.extract_variants(assembly_file, output_dir)
        self.gu.filter_SNPs(assembly_file, "filtered_snps.vcf", output_dir, params)
        self.gu.filter_Indels(assembly_file, "filtered_indels.vcf", output_dir, params)
        self.gu.exclude_filtered_variants(output_dir)
        self.gu.base_quality_score_recalibration(assembly_file, "recal_data.table", output_dir)
        self.gu.apply_BQSR(assembly_file, "recal_data.table", output_dir)
        self.gu.base_quality_score_recalibration(assembly_file, "post_recal_data.table", output_dir)
        self.gu.apply_BQSR(assembly_file,  "post_recal_data.table", output_dir)
        self.gu.filter_SNPs(assembly_file, "filtered_snps_final.vcf", output_dir, params)
      
        #Todo: To save indels also using VariationUtils or merge with snps and sort them with chr & pos and save using variaiotiontuils.
        #Todo: To get an example for saving structural variants(specially CNV) and compare with standard vcf output.

        self.gu.filter_Indels(assembly_file, "filtered_indels_final.vcf", output_dir, params)

        '''
        os.system("grep   '##fileformat' " + output_dir + "/filtered_snps_final.vcf > " + output_dir + "/sample.vcf")
        cmd = "grep -v  '##' " + output_dir + "/filtered_snps_final.vcf >> " + output_dir + "/sample.vcf"
        os.system(cmd)            # TODO : need to remove system command after fixing variationUtils.
        '''

        #Todo : check existence of final filtered finals snps.
        #Todo : chnage assembly_or_genome_ref to genome_or_assembly_ref

        #Todo: to derive name of sample_attribute_name from sample set ref by prefixing/suffixing. Attribute mapping should have one sample.
  
        save_variation_params = {'workspace_name': params['workspace_name'],
            'genome_or_assembly_ref': params['assembly_or_genome_ref'],      
            'sample_set_ref':params['input_sample_set'],
            'sample_attribute_name':'sample_attr',
            'vcf_staging_file_path': output_dir + "/filtered_snps_final.vcf",
            'variation_object_name': params['variation_object_name']
            } 

        self.vu.save_variation_from_vcf(save_variation_params)

        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created': [],
                                                'text_message': 'Success'},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_kb_GATK

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_GATK return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
