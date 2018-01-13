import sys
import json

import os
import requests

# ===================================================================================
from PythonServer.Primer import Primer
from PythonServer.Primer_set import Primer_set
from PythonServer.EnsemblRestClient import EnsemblRestClient

config = json.load(open('config.json'))


def transcript_data(species, symbol):
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol)
    transcripts = variants['Transcript']
    if transcripts:
        for t in transcripts:
            if t['display_name'] == symbol + "-201":
                #print '{display_name}: ==> {id}'.format(**t)
                juncArr = junctions(t)
                cDna = get_cdna(t['id'])
                cDna = cDna.replace("\n", "")
                #print cDna
                id_count = 1
                # forward
                optional_primers = get_optional_primers(cDna, juncArr, id_count)
                forward_primers = primer_tests(optional_primers)
                #for primer in forward_primers:
                 #   primer.printPrimer()

                # reverse
                reverse_optional_primers = get_reverse_primers(cDna, forward_primers, id_count)
                reverse_primers = primer_tests(reverse_optional_primers)
                #for primer in reverse_primers:
                 #   primer.printPrimer()

                # sets
                primers_optinal_sets = get_optional_sets(forward_primers, reverse_primers)
                primers_sets = sets_tests(primers_optinal_sets)
                sets_tests(primers_sets)
                primers_sets.sort(key=get_set_score, reverse=True)
                #print len(primers_sets)
                #for set in primers_sets:
                 #   set.set_print()
                export_to_file(primers_sets,species,symbol)


# ====================== return all junctions index ========================================

def junctions(transcript):
    junctions_arr = []
    exons_len = []
    exons = transcript['Exon']
    if exons:
        for e in exons:
            #print 'start : {start} ==> end : {end}'.format(**e)
            # if e == exons[0]:  need to check again!!!
            #   exons_len.append(transcript['end']-e['start'])
            # elif e == exons[len(exons)-1]:
            #    exons_len.append(e['end'] - transcript['start'])
            # else:
            exons_len.append(e['end'] - e['start'])
    #print exons_len
    if exons_len:
        sum = 0
        for i in exons_len:
            junctions_arr.append(sum + i)
            sum = sum + i

    return junctions_arr


# ======================= get cDna of tarscript ===========================================

def get_cdna(transcriptId):
    server = "http://rest.ensembl.org"
    ext = "/sequence/id/" + transcriptId + "?type=cdna"
    r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.text


# ======================= get all optional primers array ===========================================

def get_optional_primers(cdna, junctionArray, id_count):
    len_range = [14, 15, 16, 17, 18, 19, 20]
    threshold = [0.4, 0.5, 0.6]
    optional_primers = []
    for th in threshold:
        for l in len_range:
            for index in junctionArray:
                if ((index - (l * th)) >= 0) and ((index + (l * (1 - th))) < len(cdna)):
                    f = int(round((index - (l * th)), 0))
                    t = int(round((index + (l * (1 - th))), 0))
                    tmp_primer = Primer(id_count, "forward", cdna[f:t], f)
                    id_count += 1
                    optional_primers.append(tmp_primer)
                if th != 0.5:
                    if ((index + (l * (1 - th))) >= 0) and ((index - (l * th)) < len(cdna)):
                        f = int(round((index - (l * (1 - th))), 0))
                        t = int(round((index + (l * th)), 0))
                        tmp_primer = Primer(id_count, "forward", cdna[f:t], f)
                        id_count += 1
                        optional_primers.append(tmp_primer)

    return optional_primers


# ======================= primers test methods ===========================================

def primer_tests(optional_primers):
    primers = []
    for primer in optional_primers:
        if tm_test(primer):
            if gc_test(primer):
                if syntax_tests(primer):
                    primers.append(primer)

    #print len(primers)
    return primers


# ======================= tm test ===========================================

def tm_test(primer):
    if (primer.primer_tm() < int(config["Tm"]["Min"])) or (primer.primer_tm() > int(config["Tm"]["Max"])):
        return False
    return True


# ======================= %GC test ===========================================

def gc_test(primer):
    if (primer.precent_gc() < int(config["GC Percent"]["Min"])) \
            or (primer.precent_gc() > int(config["GC Percent"]["Max"])):
        return False
    return True


# ======================= %syntaxTests test ===========================================

def syntax_tests(primer):
    if primer.sequence[0] == 'G':
        return False

    if "GGGG" in primer.sequence:
        return False

    if "CCCC" in primer.sequence:
        return False

    if "AAAAAA" in primer.sequence:
        return False

    consecutive_c = 0
    for index in range(len(primer.sequence)):
        if primer.sequence[index] == 'C':
            consecutive_c += 1
        else:
            consecutive_c = 0
        if consecutive_c == 2:
            return False  # primer cant have 2 consecutive C in the middle
    return True


# ======================= Primer score ============================================

def primers_score(primer_sets):
    for primer_set in primer_sets:
        primer_set.get_primers_set_score()

# ======================reverse nucleotide=====================

def reverse_nucleotide(sequence):
    tmp_str = list(sequence)
    for index in range(len(sequence)):
        if sequence[index] == 'A':
            tmp_str[index] = 'T'
        elif sequence[index] == 'T':
            tmp_str[index] = 'A'
        elif sequence[index] == 'G':
            tmp_str[index] = 'C'
        else:
            tmp_str[index] = 'G'
        sequence = "".join(tmp_str)


# ======================find reverse primers=====================

def get_reverse_primers(cdna, forward_primers, id_count):
    optional_reverse_primers = []
    cdna_length = len(cdna)
    primers_length = [14, 15, 16, 17, 18, 19, 20]
    for primer in forward_primers:
        for amplicon_length in range(50, 150):
            for primer_length in primers_length:
                if cdna_length - primer.start_index >= amplicon_length + primer_length:
                    f = primer.start_index + amplicon_length
                    t = primer.start_index + amplicon_length + primer_length
                    sequence = cdna[f:t]
                    reverse_nucleotide(sequence)
                    tmp_primer = Primer(id_count, "reverse", sequence, primer.start_index + amplicon_length,
                                        primer.id)
                    optional_reverse_primers.append(tmp_primer)
                    id_count += 1
    return optional_reverse_primers


# ====================== Create sets of primers ==================================

def get_optional_sets(forward_primers, reverse_primers):
    primer_sets = []
    for rev_primer in reverse_primers:
        for for_primer in forward_primers:
            if rev_primer.pair_id == for_primer.id:
                tmp_set = Primer_set(for_primer, rev_primer)
                primer_sets.append(tmp_set)
                break
    return primer_sets


# ====================== Sets tests ==============================================

def sets_tests(primer_optional_sets):
    primer_sets = []
    for set in primer_optional_sets:
        if set.tm_dif() <= 1:
            primer_sets.append(set)
    return primer_sets


# ==================================export to file===================================

def export_to_file(primer_sets,species,symbol):
    new_file = open("../output/primer_list_"+species+"_"+symbol+".txt", "w")
    new_file.write("PRIMERS SETS of "+symbol+"("+species+" Gene):\n")
    new_file.close()
    print "Proceed to the output folder to view the results."
    count = 0;
    for index in range(0, 100):
        primer_sets[index].write_to_file("../output/primer_list_"+species+"_"+symbol+".txt")

# ================================== return primers set score ===================================

def get_set_score(primer_set):
    return primer_set.get_primers_set_score()

# ======================= MAIN ====================================================

if __name__ == '__main__':
    # if len(sys.argv) == 3:
    #   species, symbol = sys.argv[1:]
    # else:
    #    species, symbol = 'human', 'BRAF'
    species = raw_input("Enter specie:")
    symbol = raw_input("Enter symbol:")
    transcript_data(species, symbol)


