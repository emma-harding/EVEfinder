import re
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

animal = "yellow"
prefix = "af" #MUST be two letters only

def assignIDs(csv, metadata, output):
    # importing the csv files as a pandas databases
    df = pd.read_csv(csv, sep='\t')
    df.columns = ['Accession', 'Scaffold', 'Description', 'Identity', 'Length', 'E value', 'Scaffold start', 'Scaffold end', 'Query start', 'Query end', 'Query sequence', 'queryFrame']
    df['Accession'] = df['Accession'].str.split(".",expand=True)[0]

    # converting reverse frame hits into forward direction
    df['Scaffold start'], df['Scaffold end'] = np.where(df['Scaffold start'] < df['Scaffold end'], [df['Scaffold start'], df['Scaffold end']], [df['Scaffold end'], df['Scaffold start']])

    # sorting entries by scaffold and then by start value
    df.sort_values(['Scaffold', 'Scaffold start'], inplace = True, ascending = True)

    #Looping through BLAST results to associate each with an ERV
    scaffolds = []  # list containing scaffolds that have been processed
    hitID = 0
    max = 0
    ERVid = []
    for index, row in df.iterrows():
        scaffold = row['Scaffold']
        # checking if this is a new scaffold
        if sum(scaffold in s for s in scaffolds) == 0: # if this is a new scaffold
            # begin new ERV
            hitID += 1
            scaffolds.append(row['Scaffold'])
            max = row['Scaffold end']
        else:
            # checking if the hit is overlapping our current region or within 1000NT of the start or end
            if (row['Scaffold start'] <= (max + 1000)):
                if row['Scaffold end'] > max:
                    max = row['Scaffold end']
            else:
                hitID += 1
                max = row['Scaffold end']
        ERVid.append(hitID)


    # reading in table of metadata for accession numbers
    metadf = pd.read_csv(metadata)
    # merging the BLAST results with the metadata
    df['ERVid'] = ERVid
    dfmerged = pd.merge(df, metadf, on='Accession', how='inner')
    dfmerged.sort_values(['ERVid'], inplace = True, ascending = True)
    dfmerged.to_csv(output, index=False, sep='\t')

def makeNrlistRetro(csv, output):
    dfmerged = pd.read_csv(csv, sep='\t')
    #creating the nr list with metadata
    nrList = []
    ERVnumbers = list(dfmerged['ERVid'].unique())

    for i in ERVnumbers:
        #creating dataframe for each ERV
        subdf = dfmerged[dfmerged['ERVid'] == i]
        scaffold = subdf.Scaffold.iloc[0]
        hitStart = subdf['Scaffold start'].min()
        hitEnd = subdf['Scaffold end'].max()
        generaMode = subdf['Genus'].mode().values
        avident = subdf['Identity'].mean()

        # #finding the start and stop of each retroelement gene
        gagStart = (subdf[subdf.Gene == "Gag"])['Scaffold start'].min()
        gagEnd = (subdf[subdf.Gene == "Gag"])['Scaffold end'].max()
        gagGenus = (subdf[subdf.Gene == "Gag"])['Genus'].mode().values
        polStart = (subdf[subdf.Gene == "Polymerase"])['Scaffold start'].min()
        polEnd = (subdf[subdf.Gene == "Polymerase"])['Scaffold end'].max()
        polGenus = (subdf[subdf.Gene == "Polymerase"])['Genus'].mode().values
        envStart = (subdf[subdf.Gene == "Envelope"])['Scaffold start'].min()
        envEnd = (subdf[subdf.Gene == "Envelope"])['Scaffold end'].max()
        envGenus = (subdf[subdf.Gene == "Envelope"])['Genus'].mode().values
        gag = str(gagStart > 1)
        pol = str(polStart > 1)
        env = str(envStart > 1)

        #detecting recombinants
        if (pol == "True" and env == "True"):
            recombinant = ((", ").join(envGenus) != (", ").join(polGenus))
        elif (gag == "True" and env == "True"):
            recombinant = ((", ").join(envGenus) != (", ").join(gagGenus))
        #elif (gag == "True" and pol == "True"):
        #    recombinant = ((", ").join(gagGenus) != (", ").join(polGenus))
        else:
            recombinant = False
        recombinant = str(recombinant)

        #metaList = [i, scaffold, hitStart, hitEnd, avident, (", ").join(generaMode), (", ").join(genera)]
        metaList = [i, scaffold, hitStart, hitEnd, avident, recombinant, (", ").join(generaMode), gag, gagStart, gagEnd, (", ").join(gagGenus), pol, polStart, polEnd, (", ").join(polGenus), env, envStart, envEnd, (", ").join(envGenus)]
        nrList.append(metaList)

    nrdf = pd.DataFrame(nrList)
    #nrdf.columns = ['ERVid', 'Scaffold', 'hitStart', 'hitEnd', 'identity', 'GeneraMode', 'Genera']
    nrdf.columns = ['ERVid', 'Scaffold', 'hitStart', 'hitEnd', 'identity', 'Recombinant', 'GeneraMode', 'Gag', 'gagStart', 'gagEnd', 'gagGenus', 'Pol', 'polStart', 'polEnd', 'polGenus', 'Env', 'envStart', 'envEnd', 'envGenus']
    nrdf['length'] = nrdf['hitEnd'] - nrdf['hitStart']
    nrdf.to_csv(output, index=False, sep='\t')

def makeGTF(csv, GTF):
    df = pd.read_csv(csv, sep='\t')
    gtfFrame = pd.DataFrame().assign(seqname=df['Scaffold'], source="EVEfinder", feature=df['GeneraMode'], start=df['hitStart'], end=df['hitEnd'], score=".", strand="+", frame="0", attribute=df['ERVid'])
    gtfFrame['attribute'] = gtfFrame['attribute'].astype(str)
    gtfFrame.to_csv(GTF, index=False, sep='\t')

def generateStats(csv, statsOutput):
    nrdf = pd.read_csv(csv, sep='\t')
    nrdf = nrdf.astype({'Gag': str, 'Pol': str, 'Env': str, 'Recombinant':str}, errors='raise')
    total = nrdf.ERVid.max()
    retrotransposonCount = nrdf.GeneraMode.str.count("Retrotransposon").sum()
    gagCount = nrdf.Gag.str.count("True").sum()
    polCount = nrdf.Pol.str.count("True").sum()
    envCount = nrdf.Env.str.count("True").sum()
    recombinantCount = nrdf.Recombinant.str.count("True").sum()
    betaCount = nrdf.GeneraMode.str.count("Betaretrovirus").sum()
    gammaCount = nrdf.GeneraMode.str.count("Gammaretrovirus").sum()
    alphaCount = nrdf.GeneraMode.str.count("Alpharetrovirus").sum()
    lentiCount = nrdf.GeneraMode.str.count("Lentivirus").sum()
    deltaCount = nrdf.GeneraMode.str.count("Deltaretrovirus").sum()
    epsilonCount = nrdf.GeneraMode.str.count("Epsilonretrovirus").sum()
    spumaCount = nrdf.GeneraMode.str.count("Spumavirus").sum()
    totalERVs = betaCount + gammaCount + alphaCount + lentiCount + deltaCount + epsilonCount + spumaCount
    class1 = epsilonCount + gammaCount
    class2 = betaCount + lentiCount + alphaCount + deltaCount
    class3 = spumaCount
    longest = nrdf.length.max()

    #writing the overall stats file
    with open(statsOutput, 'w') as f:
        f.write('Total: ' + str(total) + '\n' \
                + 'Retrotransposons: ' + str(retrotransposonCount) + ' ' + str((retrotransposonCount / total) * 100) + '\n' \
                + 'Endogenous retroviruses: ' + str(totalERVs) + ' ' + str((totalERVs / total) * 100) + '\n' \
                + 'Gag: ' + str(gagCount) + ' ' + str((gagCount / totalERVs) * 100) + '\n' \
                + 'Pol: ' + str(polCount) + ' ' + str((polCount / totalERVs) * 100) + '\n' \
                + 'Env: ' + str(envCount) + ' ' + str((envCount / totalERVs) * 100) + '\n' \
                + 'Recombinants: ' + str(recombinantCount) + ' ' +  str((recombinantCount / totalERVs) * 100) + '\n' \
                + 'Class 1: ' + str(class1) + '\n' \
                + 'Class 2: ' + str(class2) + '\n' \
                + 'Class 3: ' + str(class3) + '\n' \
                + 'Betaretrovirus: ' + str(betaCount) + ' ' + str((betaCount / totalERVs) * 100) + '\n' \
                + 'Gammaretrovirus: ' + str(gammaCount) + ' ' + str((gammaCount / totalERVs) * 100) + '\n' \
                + 'Alpharetrovirus: ' + str(alphaCount) + ' ' + str((alphaCount / totalERVs) * 100) + '\n' \
                + 'Lentivirus: ' + str(lentiCount) + ' ' + str((lentiCount / totalERVs) * 100) + '\n' \
                + 'Deltavirus: ' + str(deltaCount) + ' ' + str((deltaCount / totalERVs) * 100) + '\n' \
                + 'Epsilonretrovirus: ' + str(epsilonCount) + ' ' + str((epsilonCount / totalERVs) * 100) + '\n' \
                + 'Spumavirus: ' + str(spumaCount) + ' ' + str((spumaCount / totalERVs) * 100) + '\n' \
                + 'Longest: ' + str(longest))
    f.close()

def extractERVs(csv, offset, genome, fasta, prefix):
    print("indexing genome")
    genome = SeqIO.index(genome, "fasta")
    ERVs = []
    oldScaffold = []
    df = pd.read_csv(csv, sep='\t')
    for index, row in df.iterrows():
        print(row['ERVid'])
        if("|" in row['Scaffold']):
            scaffold = (row['Scaffold']).split('|')[1]
        else:
            scaffold = row['Scaffold']
        start = int(row['hitStart'] - offset)
        stop = int(row['hitEnd'] + offset)
        if(scaffold != oldScaffold):
            scaffoldSeq = genome[scaffold].seq
            oldScaffold = scaffold
        currentERV = SeqRecord(scaffoldSeq[start-1:stop])
        currentERV.id = prefix + str(row['ERVid'])
        currentERV.description = str(scaffold) + ":" + str(start) + "-" + str(stop)
        #SeqIO.write(currentERV, fasta+' '+str(offset)+'_ERVs.fasta', "fasta")
        ERVs.append(currentERV)
        SeqIO.write(ERVs, (fasta+'_'+str(offset)+'_ERVs.fasta'), "fasta")

# CHECK/CHANGE THESE BEFORE EACH RUN OF THIS SCRIPT
BlastOutput = "BLAST/" + animal + '_tblastn.tsv' # output of initial genome blast
ClassifiedBlastOutput = "BLAST/" + animal + '_tblastn_id.tsv' # non-redundant list of scaffold regions with hits
nrListOutput = "nrList/" + animal + '_nrList.tsv'
genome = "D:\PhD\Viral_discovery\Samples\Genomes\GCA_016432865.1_AdamAnt_genomic.fna" # genome
metadata = "RefRetro_metadata.csv"
statFile = "nrList/" + animal + '_stats.txt'
fasta = "ERVs/" + animal
Longfasta = "ERVs/" + animal + '_longERVs.fasta'
genBank = "ERVs/" + animal
GTF = animal + '.gtf'

# DO NOT CHANGE THINGS BELOW THIS LINE
if len(prefix) == 2:
    assignIDs(BlastOutput, metadata, ClassifiedBlastOutput)
    makeNrlistRetro(ClassifiedBlastOutput, nrListOutput)
    makeGTF(nrListOutput, GTF)
    generateStats(nrListOutput, statFile)
    extractERVs(nrListOutput, 0, genome, fasta, prefix)

else:
    print('Prefix is incorrect length.')
    exit()


