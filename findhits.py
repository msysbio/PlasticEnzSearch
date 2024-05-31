from Bio.Blast import NCBIXML

def get_best_hits(blast_output_file, num_hits=3):
    with open(blast_output_file, 'r') as file:
        blast_records = NCBIXML.parse(file)

        best_hits = []
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    best_hits.append((alignment, hsp))

        best_hits.sort(key=lambda x: x[1].bits, reverse=True)
        best_hits = best_hits[:num_hits]

        return best_hits

def print_hits(xmlfile):
    best_hits = get_best_hits(xmlfile)
    if best_hits:
        for i, (alignment, hsp) in enumerate(best_hits, start=1):
            print(f"Hit {i}: {alignment.title}")
        for i, (alignment, hsp) in enumerate(best_hits, start=1):
            print(f"Hit {i} & {alignment.length} & {hsp.bits} & {hsp.score} & {hsp.expect} & {hsp.identities}/{hsp.align_length} & {hsp.gaps}/{hsp.align_length} & {hsp.query_start}-{hsp.query_end} & {hsp.sbjct_start}-{hsp.sbjct_end} \\\\ \\hline")
    else:
        print("No hits found.")

print_hits('/home/jasper/Thesis/PET-EXPOSED/output/Daphnia_PET_pbat_hmm_output.fasta.xml')


