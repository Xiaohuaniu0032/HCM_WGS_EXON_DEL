perl extract_mane_exon_bed.pl \
  --gtf hg19.ncbiRefSeq.gtf \
  --mane ../canonical_transcripts/RefSeq_MANE_Select.xls \
  --out RefSeq_MANE_Select.exon.bed \
  --unmatched RefSeq_MANE_Select.unmatched.tsv \
  2> log
