import logging
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

def blast_local(query, output="/home/jasper/Thesis/blast/blast_results.xml"):
    # path to local blast binary
    blastp_path = "/home/jasper/Thesis/blast/ncbi-blast-2.15.0+/bin/blastp"

    # path to local blast database
    database_path = "/home/jasper/Thesis/blast/database/swissprot"

    # create a blast query command
    blastp_cline = NcbiblastpCommandline(cmd=blastp_path, query=query, db=database_path, evalue=0.01, outfmt=5, out=output)

    # execute the command
    stdout, stderr = blastp_cline()

    # parse the XML output
    result_handle = open(output)
    
    return result_handle

def blast_parse(record, local=True):
    try:
        if local:
            result_handle = blast_local(record.seq) 
        else:
            result_handle = NCBIWWW.qblast("blastp", "swissprot", record.seq, expect=0.01, filter=True, word_size=6)

        #parse the result
        blast_record = NCBIXML.read(result_handle)
        return blast_record
    
    except Exception as e:
        logging.error(f"Error: {e}")
