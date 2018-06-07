from __future__ import division
import sys
import json
import requests

# ===================================================================================
from Primer import Primer
from Primer_set import Primer_set
from EnsemblRestClient import EnsemblRestClient

config = json.load(open('config.json'))





def transcript_data(species, symbol):
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol)
    transcripts = variants['Transcript']
    if transcripts:
        for t in transcripts:
            if t['display_name'] == symbol + "-201":
                juncArr = junctions(t)
                cDna = get_cdna(t['id'])
                cDna = cDna.replace("\n", "")
                cDna = cDna.replace(">"+t['id'], "")
                cDna = remove_utr(cDna, t['UTR'])
                utr_len = get_start_utr_length(t['UTR'])
                juncArr = fix_juncArr_after_UTR_remove(juncArr, utr_len)
                id_count = 1
                id_count_set = 1
                print(cDna)
                # On junction options
                print "Searching for FORWARD Primers..."
                junc_optional_forward = get_optional_primers(cDna, juncArr, id_count, "forward")


                # reverse
                print "Searching for REVERSE Primers..."
                # on junction reverse
                junc_optional_reverse = get_optional_primers(cDna, juncArr, id_count, "reverse")

                cr = 0
                for rprimer in junc_optional_reverse:
                    if rprimer.id == 4163:
                        cr += 1
                cf = 0
                for fprimer in junc_optional_forward:
                    if fprimer.id == 4126:
                        cf += 1

                # other
                #reverse_primers = get_reverse_primers(cDna, junc_optional_forward, id_count)

                # sets
                print "Creating Primers pairs and giving scores..."
                on_junk_sets = junktion_sets(junc_optional_forward, junc_optional_reverse, id_count_set)
                on_junk_sets.sort(key=get_set_score, reverse=True)
                #primers_sets = get_optional_sets(junc_optional_forward, reverse_primers)
                #primers_sets.sort(key=get_set_score, reverse=True)


                print "Export to output file..."
                #export_to_file(on_junk_sets,primers_sets,species,symbol)
                export_to_file(on_junk_sets, species, symbol)


# ====================== remove UTR ========================================


def remove_utr(cdna,utr_object):
    start_utr = (utr_object[0]['end']-utr_object[0]['start'])+1#include the letter in the end
    end_utr = (utr_object[1]['end'] - utr_object[1]['start']) + 1  # include the letter in the end
    #print start_utr
    #print end_utr
    return cdna[start_utr:len(cdna)-end_utr]
# ====================== start UTR length ========================================


def get_start_utr_length(utr_object):
    start_utr = (utr_object[0]['end'] - utr_object[0]['start']) + 1  # include the letter in the end
    return start_utr
# ====================== fix juncArr ========================================


def fix_juncArr_after_UTR_remove(juncArr, utr_len):
    fix_juncArr = []
    for i in juncArr:
        if i-utr_len <= 0:
            continue
        fix_juncArr.append(i-utr_len)
    return fix_juncArr
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
            exons_len.append(e['end'] - e['start']+1)
    #print exons_len
    cdnalen = 0
    for i in exons_len:
        cdnalen = cdnalen+i
    cdnalen = cdnalen + len(exons_len)
    #print cdnalen
    if exons_len:
        sum = 0
        for i in exons_len:
            junctions_arr.append(sum + i)
            sum = sum + i
    #print junctions_arr
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

def get_optional_primers(cdna, junctionArray, id_count ,kind):
    config = json.load(open("config.json"))
    optional_primers = []   # list of primers
    # get the sequence length range from config file
    len_range = []
    for i in range(int(config["Length"]["Min"]), int(config["Length"]["Max"])+1):
        len_range.append(i)

    for l in len_range:
        # get the junction percentage range from config file
        min_threshold = int(round(l * float(config["Threshold"]["Min"]), 0))
        max_threshold = int(round(l * float(config["Threshold"]["Max"]), 0))
        threshold = range(min_threshold, max_threshold + 1)
        for th in threshold:
            for index in junctionArray:
                if ((index - th) >= 0) and ((index + (l-th)) < len(cdna)):
                    f = index - th
                    t = f+l
                    if kind == "forward":
                        tmp_primer = Primer(id_count, "forward", cdna[f:t], f, round(th/l, 2))
                    else:
                        rev_seq = cdna[f:t]
                        rev_seq = reverse_nucleotide(rev_seq)
                        tmp_primer = Primer(id_count, "reverse", rev_seq, f, round(th/l, 2))
                    id_count += 1
                    optional_primers.append(tmp_primer)
                    '''
                if th != l*0.5:       # add both sides , to avoid duplicate primer on half len
                    if ((index - (l - th)) >= 0) and ((index + th) < len(cdna)):
                        f = index - (l-th)
                        t = f + l
                        if cdna[f:t] == reverse_nucleotide("TGAGGTGTAGGTGCTGTCACA"):
                            pass
                        if kind == "forward":
                            tmp_primer = Primer(id_count, "forward", cdna[f:t],  f, round(th/l, 2))
                        else:
                            rev_seq = cdna[f:t]
                            rev_seq = reverse_nucleotide(rev_seq)
                            tmp_primer = Primer(id_count, "reverse", rev_seq,  f, round(th/l, 2))
                        id_count += 1
                        optional_primers.append(tmp_primer)'''

    optional_primers = primer_tests(optional_primers)
    return optional_primers


# ======================= primers test methods ===========================================

def primer_tests(optional_primers):
    primers = []
    for primer in optional_primers:
        if primer.sequence == "TGGACCAGGAGGAAGCTCAG":
            pass
        if tm_test(primer):
            if gc_test(primer):
                #if syntax_tests(primer):
                if palindrome_test(primer):
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

# ======================= palindrome test ===========================================

def palindrome_test(primer):
    if abs(primer.palindrome_length) > 10:
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
    sequence = sequence[::-1]
    return sequence




# ======================find reverse primers=====================

def get_reverse_primers(cdna, forward_primers, id_count):
    optional_reverse_primers = []
    cdna_length = len(cdna)
    primers_length = []
    amplicon_length = []
    for i in range(int(config["Length"]["Min"]), int(config["Length"]["Max"])+1):
        primers_length.append(i)
    for i in range(int(config["Amplicon Length"]["Min"]), int(config["Amplicon Length"]["Max"])+1):
        amplicon_length.append(i)
    for primer in forward_primers:
        for amplic_len in amplicon_length:
            for primer_length in primers_length:
                if cdna_length - primer.start_index >= amplic_len + primer_length:
                    f = primer.start_index + amplic_len
                    t = primer.start_index + amplic_len + primer_length
                    sequence = cdna[f:t]
                    sequence = reverse_nucleotide(sequence)
                    tmp_primer = Primer(id_count, "reverse", sequence, primer.start_index + amplic_len, pair_id=primer.id)
                    optional_reverse_primers.append(tmp_primer)
                    id_count += 1

    optional_reverse_primers = primer_tests(optional_reverse_primers)
    return optional_reverse_primers


# ====================== Create sets of primers ==================================

def get_optional_sets(forward_primers, reverse_primers, id_count_set):
    primer_sets = []
    for rev_primer in reverse_primers:
        for for_primer in forward_primers:
            if rev_primer.pair_id == for_primer.id:
                tmp_set = Primer_set(for_primer, rev_primer, id_count_set)
                id_count_set += 1
                primer_sets.append(tmp_set)
                break

    #primer_sets = sets_tests(primer_sets)       # send to tests
    return primer_sets


# ====================== Sets tests ==============================================

def sets_tests(primer_optional_sets):
    tm_diff_range = []
    #getting tm difference ranges from config file
    tm_dif_min = int(config["Temperature difference"]["Min"])
    tm_dif_max = int(config["Temperature difference"]["Max"])
    primer_sets = []
    for set in primer_optional_sets:
        if set.forward_primer.sequence == "TGGACCAGGAGGAAGCTCAG" and set.reverse_primer.sequence == "TGTCCGCTGTGCCAGATATG":
            pass
        if tm_dif_min <= set.tm_dif() <= tm_dif_max:
            primer_sets.append(set)
    return primer_sets


# ============== Check for sets of forward&reverse on junktion  ==============================================

def junktion_sets(forward_on_junk, reverse_on_junk, id_count_set):
    on_junc_sets = []
    for forward_primer in forward_on_junk:
        if forward_primer.sequence == "TGGACCAGGAGGAAGCTCAG":
            pass
        for reverse_primer in reverse_on_junk:
            if reverse_primer.sequence == "TGTCCGCTGTGCCAGATATG":
                pass
            if (reverse_primer.start_index+reverse_primer.length)-forward_primer.start_index in range((int(config["Amplicon Length"]["Min"])), (int(config["Amplicon Length"]["Max"]))+1):
                set = Primer_set(forward_primer, reverse_primer, id_count_set)
                id_count_set += 1
                on_junc_sets.append(set)

    on_junc_sets = sets_tests(on_junc_sets)
    return on_junc_sets






# ==================================export to file===================================

def export_to_file(on_junk_sets,species,symbol):
    #new_file_junk = open("../output/primer_list_junktion_"+species+"_"+symbol+".txt", "w")
   # new_file_junk.write("Forward & Reverse on junktion PRIMERS SETS of "+symbol+"("+species+" Gene):\n")
    #new_file_junk.close()

    new_file = open("../output/primer_list_" + species + "_" + symbol + ".txt", "w")
    new_file.write("\n"+
                   "*************************************************************************************************\n"+
                   "********************Forward & Reverse on junktion PRIMERS SETS of "+symbol+"("+species+" Gene):**************\n"+
                   "*************************************************************************************************\n\n")
    new_file.close()

    max_range = 1000
    if len(on_junk_sets) < max_range:
        max_range = len(on_junk_sets)
    for index in range(0, max_range):        # On junk sets
        if on_junk_sets[index].id == 49796:
            pass
        on_junk_sets[index].write_to_file("../output/primer_list_"+species+"_"+symbol+".txt")

    '''max_range = 15
    if len(primer_sets) < max_range:
        max_range = len(primer_sets)
    with open("../output/primer_list_" + species + "_" + symbol + ".txt", "a") as f:
        f.write("\n"+
        "******************************************************************************************************\n"+
                "********************************PRIMERS SETS of " + symbol + "(" + species + " Gene):*************************************\n"+
                "******************************************************************************************************\n\n")
    for index in range(0, max_range):
        primer_sets[index].write_to_file("../output/primer_list_"+species+"_"+symbol+".txt")
    '''
    print "Done! \nProceed to the output folder to view the results."


# ================================== return primers set score ===================================

def get_set_score(primer_set):
    return primer_set.get_primers_set_score()
# ================ method to run between 2 float numbers with float jumping ======================
def frange(x, y, jump):
  while x < y:
    yield x
    x += jump
