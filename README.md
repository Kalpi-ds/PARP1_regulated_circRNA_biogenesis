# PARP1_regulated_circRNA_biogenesis_2023

PARP1 Regulates Circular RNA Biogenesis though Control of Transcriptional Dynamics

Eleazer, Rebekah, Kalpani De Silva, Kalina Andreeva, Zoe Jenkins, Nour Osmani, Eric C. Rouchka, and Yvonne Fondufe-Mittendorf. 2023. "PARP1 Regulates Circular RNA Biogenesis though Control of Transcriptional Dynamics" Cells 12, no. 8: 1160. https://doi.org/10.3390/cells12081160

circRNA detection

miRNA analysis

1. once seekCRIT done. Divid circRNAs across samples --->parseOutput.pl
2. get the fasta sequences for the divided circRNAs ---->parseSeekCRIT_2.pl
3. find miRNA seeds in these fastas --->parseMIRsV2.pl
MIRBASE_22.1_hsa_mature.fa get this from mibase it contains mature miRNA seeds for human.
4. find the enriched, most frequent enriched miRNAs --->process_MostFrequentenrichedMIRs_list.R
5. get EntrezIDs for the target gens 0f these miRNAs --->extract_EntrezIDs_for_miRNATarget_Genes.R
hsa_MTI.txt get this from miRTarbase. it is a list of miRNa target genes for human
6. Use these Entrez for GO-clusterprofiler, KEGG-grofiler --->ClusterProfiler.R, gprofiler.R, enrichmentBarPlot.R 
7. C:\Kalpani\KY-INBRE PostDoc\Yvonne\Bekah\2024\compareiAsOVESaTBshSATB2onlymiRNAs Compare miRNas between conditions --> experimentalandcontrolOnlymiRNAOverlap.R
8. C:\Kalpani\KY-INBRE PostDoc\Yvonne\Bekah\2024\compareiAsOVESaTBshSATB2onlycircRnAs ompare circRNAs between conditiond ---> experimentalConditionOnlycircRNAOverlap.R
