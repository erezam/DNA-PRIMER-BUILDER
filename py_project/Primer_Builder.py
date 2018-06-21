
import sys
import json
import requests

# ===================================================================================
from Primer import Primer
from Primer_set import Primer_set
from EnsemblRestClient import EnsemblRestClient
import TmPredictor
import ScorePredictor
# load config file with the parameters ranges
config = json.load(open('config.json'))


# ====================== Build ========================================
# main func of Primer_Builder that manage the full process
def build(species, symbol):
    # init tm predictor once to reset the global variable regr for tm calc with logistic reg
    ScorePredictor.init()
    TmPredictor.init()
    # get the needed "201" transcript
    transcript_201 = transcript_data(species, symbol)
    juncArr = junctions(transcript_201)
    cDna = get_cdna(transcript_201['id'])
    cDna = cDna.replace("\n", "")
    cDna = cDna.replace(">" + transcript_201['id'], "")
    cDna = remove_utr(cDna, transcript_201['UTR'])
    # UTR CUT, utr is a unnecessary range of the cdna
    utr_len = get_start_utr_length(transcript_201['UTR'])
    juncArr = fix_juncArr_after_UTR_remove(juncArr, utr_len)

    # id counters for primers/primers_sets
    id_count = 1
    id_count_set = 1

    # On junction primers !
    # forward
    print("Searching for FORWARD Primers...")
    junc_optional_forward = get_optional_primers(cDna, juncArr, id_count, "forward")

    # reverse
    print("Searching for REVERSE Primers...")
    # on junction reverse
    junc_optional_reverse = get_optional_primers(cDna, juncArr, id_count, "reverse")

    # not on junctions reverse primers
    # reverse_primers = get_reverse_primers(cDna, junc_optional_forward, id_count)

    # creating sets
    print("Creating Primers pairs and giving scores...")
    # on junction sets
    on_junk_sets = junktion_sets(junc_optional_forward, junc_optional_reverse, id_count_set)
    on_junk_sets.sort(key=get_set_score, reverse=True)
    # other sets
    # primers_sets = get_optional_sets(junc_optional_forward, reverse_primers)
    # primers_sets.sort(key=get_set_score, reverse=True)

    print("Export to output file...")
    # export_to_file(on_junk_sets,primers_sets,species,symbol)
    export_to_file(on_junk_sets, species, symbol)


# ====================== transcript data ========================================
# init the rest client and get back the full json object for the specific gene
def transcript_data(species, symbol):
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol)
    transcripts = variants['Transcript']
    if transcripts:
        for t in transcripts:
            if t['display_name'] == symbol + "-201":
                return t


# ======================= get cDna of transcript ===========================================

def get_cdna(transcript_id):
    server = "http://rest.ensembl.org"
    ext = "/sequence/id/" + transcript_id + "?type=cdna"
    r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.text

# ====================== remove UTR ========================================
# remove the utr from the cDna
def remove_utr(cdna,utr_object):
    start_utr = (utr_object[0]['end']-utr_object[0]['start'])+1#include the letter in the end
    end_utr = (utr_object[1]['end'] - utr_object[1]['start']) + 1  # include the letter in the end
    return cdna[start_utr:len(cdna)-end_utr]


# ====================== start UTR length ========================================

def get_start_utr_length(utr_object):
    start_utr = (utr_object[0]['end'] - utr_object[0]['start']) + 1  # include the letter in the end
    return start_utr


# ====================== fix juncArr ========================================
# fix the junctions index arr after utr remove
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
            exons_len.append(e['end'] - e['start']+1)
    if exons_len:
        sum = 0
        for i in exons_len:
            junctions_arr.append(sum + i)
            sum = sum + i
    return junctions_arr


# ======================= get all optional primers array on junctions ===========================================

def get_optional_primers(cdna, junctionArray, id_count, kind):
    # this func find all the primers that sits on junctions by kind(forward/reverse)
    #test
    #test_forward_primer = "TGGACCAGGAGGAAGCTCAG"
    #test_reverse_primer = "TGTCCGCTGTGCCAGATATG"

    optional_primers = []   # list of primers
    # get the sequence length range from config file
    len_range = []
    for i in range(int(config["Length"]["Min"]), int(config["Length"]["Max"])+1):
        len_range.append(i)

    for l in len_range:
        # get the junction percentage range from config file
        min_threshold = int(round(l * float(config["Threshold"]["Min"]), 0))
        max_threshold = int(round(l * float(config["Threshold"]["Max"]), 0))
        tmp_primer = None

        threshold = list(range(min_threshold, max_threshold + 1))
        for th in threshold:
            for index in junctionArray:
                if ((index - th) >= 0) and ((index + (l-th)) < len(cdna)):
                    f = index - th
                    t = f+l
                    if kind == "forward":
                        #test
                        #if cdna[f:t] == test_forward_primer:
                        tmp_primer = Primer(id_count, "forward", cdna[f:t], f, round(th/l, 2))
                    else:
                        rev_seq = cdna[f:t]
                        rev_seq = reverse_nucleotide(rev_seq)
                        #test
                        #if rev_seq == test_reverse_primer:
                        tmp_primer = Primer(id_count, "reverse", rev_seq, f, round(th/l, 2))
                    id_count += 1
                    #if tmp_primer != None:
                    optional_primers.append(tmp_primer)
                    #    tmp_primer = None


    optional_primers = primer_tests(optional_primers)
    return optional_primers


# ======================= primers test methods ===========================================
# remove all bad primers that dont fits the config parmeters
def primer_tests(optional_primers):
    primers = []
    for primer in optional_primers:
        if tm_test(primer):
            if gc_test(primer):
                #if syntax_tests(primer):
                if palindrome_test(primer):
                    primers.append(primer)
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
    if abs(primer.palindrome_length) >= 7:
        return False
    return True

# ======================= %syntaxTests test ===========================================
# not in use!! this test not fit the biologist demands
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


# ======================reverse nucleotide=====================
# turn the nucletide to rev (exmp : "GCA" to "TGC")
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
# not in use, finds rev primers not on junk by amplicon len
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
# not in use for now , create pairs of primers (Not necessarily on junc)
def get_optional_sets(forward_primers, reverse_primers, id_count_set):
    primer_sets = []
    for rev_primer in reverse_primers:
        for for_primer in forward_primers:
            if rev_primer.pair_id == for_primer.id:
                tmp_set = Primer_set(for_primer, rev_primer, id_count_set)
                id_count_set += 1
                primer_sets.append(tmp_set)
                break

    primer_sets = sets_tests(primer_sets)       # send to tests
    return primer_sets


# ====================== Sets tests ==============================================
# remove all the bad sets
def sets_tests(primer_optional_sets):
    #getting tm difference ranges from config file
    tm_dif_min = int(config["Temperature difference"]["Min"])
    tm_dif_max = int(config["Temperature difference"]["Max"])
    primer_sets = []
    for set in primer_optional_sets:
        if tm_dif_min <= set.tm_dif() <= tm_dif_max:
            primer_sets.append(set)
    return primer_sets


# ============== Check for sets of forward&reverse on junction  ==============================================
# create the sets of primers on junctions
def junktion_sets(forward_on_junk, reverse_on_junk, id_count_set):
    on_junc_sets = []
    for forward_primer in forward_on_junk:
        for reverse_primer in reverse_on_junk:
            # check amplicon len
            if (reverse_primer.start_index+reverse_primer.length)-forward_primer.start_index in range((int(config["Amplicon Length"]["Min"])), (int(config["Amplicon Length"]["Max"]))+1):
                set = Primer_set(forward_primer, reverse_primer, id_count_set)
                id_count_set += 1
                on_junc_sets.append(set)

        ''' Add primers to json data
    new_file = open("../output/data.txt", "w")
    for set in on_junc_sets:
        if set.id % 5000 == 0:
            new_file.write("{\n\t'F_seq':'%s',\n\t'F_tm':'%.2f',\n\t'F_gc':'%s',\n\t'F_len':'%s',\n\t'F_pali':'%s',\n\t"
                           "'R_seq':'%s',\n\t'R_tm':'%.2f',\n\t'R_gc':'%s',\n\t'R_len':'%s',\n\t'R_pali':'%s',\n\t"
                           "'Tm_dif':'%s',\n\t'Amp_len':'%s',\n\t'score':%s\n},\n" % (
                           set.forward_primer.sequence, set.forward_primer.primer_tm(), set.forward_primer.precent_gc(),
                           set.forward_primer.length, set.forward_primer.palindrome_length,
                           set.reverse_primer.sequence, set.reverse_primer.primer_tm(), set.reverse_primer.precent_gc(),
                           set.reverse_primer.length, set.reverse_primer.palindrome_length,
                           set.tm_dif(), set.get_amplicon_length(), 0))

    new_file.close()
    '''
    on_junc_sets = sets_tests(on_junc_sets)
    return on_junc_sets


# ================================== return primers set score ===================================

def get_set_score(primer_set):
    return primer_set.get_primers_set_score()


# ==================================export to file===================================
# write the best 100 primers to output file
def export_to_file(on_junc_sets, species, symbol):

    new_file = open("../output/primer_list_" + species + "_" + symbol + ".txt", "w")
    new_file.write("\n**********************************************************************************\
                    \n********* Forward & Reverse on junction PRIMERS SETS of "+ symbol + "(" + species+" Gene): *********\
                    \n**********************************************************************************\n\n")
    new_file.close()

    max_range = 100
    if len(on_junc_sets) < max_range:
        max_range = len(on_junc_sets)
    for index in range(0, max_range):        # On junk sets
        if on_junc_sets[index].id == 49796:
            pass
        on_junc_sets[index].write_to_file("../output/primer_list_" + species + "_" + symbol + ".txt")
    for set in on_junc_sets:
        if set.id == 311432:
            set.write_to_file("../output/primer_list_"+species+"_"+symbol+".txt")

    print("Done! \nProceed to the output folder to view the results.")


