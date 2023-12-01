# EVEfinder
Version of EVEfinder used for "Invasion and amplification of endogenous retroviruses in Dasyuridae marsupials" (Nov 2023).

Workflow to find EVEs from animal genomes

This python program (EVEfinder) takes tab-delimited BLAST output (Retrovirus AA query search against genome) and runs through five sequential functions:

    1. assignIDs: this function classifies each BLAST hit into ERVs based on proximity. Hits within 1000NT of each other are classified as a single ERV.
    2. makeNrlistRetro: this function generates a tsv file of unique ERVs with information including their location, length, retroviral classification and genes present (Gag/Pol/Env). This function also detects if retroviruses are recombinant if the Pol and Env genera are different.
    3. makeGTF: converts the nrList into a GTF file of ERV annotations
    4. generateStats: this function produces a txt file with overall statistics of the ERVs. Note it does double-count ERVs, especially when calculating the genera represented.
    5. extractERVs: this function extracts ERV sequences as a .fasta file from the genome.
    6. addAnnotations: this function converts the ERV fasta sequences into GenBank format and adds annotations for the gag/pol/env genes based on the initial BLAST hits.

