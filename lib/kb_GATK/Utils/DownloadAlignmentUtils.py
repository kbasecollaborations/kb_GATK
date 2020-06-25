from installed_clients.ReadsAlignmentUtilsClient import ReadsAlignmentUtils

class DownloadAlignmentUtils:

    def __init__(self, callback_url):
       pass 

    def downloadreadalignment(self, source_ref, params, callback_url):
        self.callback_url = callback_url
        self.ru = ReadsAlignmentUtils(self.callback_url)
        params['source_ref'] = source_ref
        params['downloadSAM'] = 1

        params['destination_dir'] = '/kb/module/work/tmp'
        params['stats'] = {
                    "properly_paired":1,
                    "multiple_alignments":1,
                    "singletons":1,
                    "alignment_rate":1,
                    "unmapped_reads":1,
                    "mapped_reads":1,
                    "total_reads":1
                 }
        self.ru.download_alignment(params)
  


