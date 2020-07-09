/*
A KBase module: kb_GATK
*/

module kb_GATK {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */

    typedef structure {
        string alignment_ref;
        string assembly_or_genome_ref;
        string output_variant_object;
        string snp_qd_filter;
        string snp_fs_filter;
        string snp_mq_filter;
        string snp_sor_filter;
        string snp_mqrankSum_filter;
        string snp_readposranksum_filter;
        string indel_qd_filter;
        string indel_fs_filter;
        string indel_mq_filter;
        string indel_sor_filter;
        string indel_mqrankSum_filter;
        string indel_readposranksum_filter;
    } Inparams;

    funcdef run_kb_GATK(Inparams params) returns (ReportResults output) authentication required;

};
