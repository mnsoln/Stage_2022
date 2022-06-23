#python3.9 pos.py --fasta /Users/smedina/Desktop/Stage/test/TAIR10_Chr.all.fasta --gff /Users/smedina/Desktop/Stage/test/Araport11_GFF3_genes_transposons.201606.gff

#python3.9 data.py --reference /Users/smedina/Documents/GitHub/Stage_2022/RefPosC.txt --bismark /Users/smedina/Desktop/Stage/test/mcseq_10_Col_0_rep2_R1_val_1_bismark_pe.deduplicated.CX_report.txt

#python3.9 pos.py --fasta /Users/smedina/Documents/BdistachyonBd21_3_537_v1.0.fa --gff /Users/smedina/Documents/BdistachyonBd21_3_537_v1.2.gene.gff3

#python3.9 data.py --reference /Volumes/DATA/GitHub/Stage_2022/ReferencePosC.txt --bismark /Volumes/DATA/Bd05_1_val_1_bismark_hisat2_pe.deduplicated.CX_report.txt --debug

python3.9 data.py --reference /Volumes/DATA/GitHub/Stage_2022/ReferencePosC.txt --bismark /Volumes/DATA/Bismark_from_Deepsignal/BD21_rep1.deepsignal_to_bismark.tsv --debug
