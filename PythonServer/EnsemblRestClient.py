from __future__ import division
import sys
import urllib
import urllib2
import json
import time
import requests

# ===================================================================================
from PythonServer.Primer import Primer


class EnsemblRestClient(object):
    # ===================== init ====================================================
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    # =================== database server connection ================================
    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urllib.urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = urllib2.Request(self.server + endpoint, headers=hdrs)
            response = urllib2.urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except urllib2.HTTPError, e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    # ===================== get variants query =============================================
    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            '/xrefs/symbol/{0}/{1}'.format(species, symbol),
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/lookup/id/{0}'.format(stable_id),
                params={'expand': '1'}
            )
            return variants
        return None

# ======================== get transcript data =============================================


def transcript_data(species, symbol):
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol)
    transcripts = variants['Transcript']
    if transcripts:
        for t in transcripts:
            if t['display_name'] == symbol+"-201":
                print '{display_name}: ==> {id}'.format(**t)
                juncArr=junctions(t)
                cDna = get_cdna(t['id'])
                cDna = cDna.replace("\n","")
                print cDna
                optional_primers = get_optional_primers(cDna, juncArr)
                primers = primer_tests(optional_primers)
                primers_score(primers)
                for primer in primers:
                    primer.printPrimer()

# ====================== return all junctions index ========================================


def junctions(transcript):
    junctions_arr = []
    exons_len = []
    exons = transcript['Exon']
    if exons:
        for e in exons:
            print 'start : {start} ==> end : {end}'.format(**e)
            #if e == exons[0]:  need to check again!!!
            #   exons_len.append(transcript['end']-e['start'])
            #elif e == exons[len(exons)-1]:
            #    exons_len.append(e['end'] - transcript['start'])
            #else:
            exons_len.append(e['end'] - e['start'])
    print exons_len
    if exons_len:
        sum = 0
        for i in exons_len:
            junctions_arr.append(sum+i)
            sum = sum+i

    return junctions_arr


# ======================= get cDna of tarscript ===========================================

def get_cdna(transcriptId):
    server = "http://rest.ensembl.org"
    ext = "/sequence/id/"+transcriptId+"?type=cdna"
    r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.text

# ======================= get all optional primers array ===========================================


def get_optional_primers(cdna,junctionArray):
    len_range=[14, 15, 16, 17, 18, 19, 20]
    threshold = [0.4, 0.5, 0.6]
    optional_primers = []
    for th in threshold:
        for l in len_range:
            for index in junctionArray:
                if ((index-(l*th))>=0) and ((index+(l*(1-th))) < len(cdna)):
                    f = int(round((index-(l*th)), 0))
                    t = int(round((index+(l*(1-th))), 0))
                    tmp_primer = Primer("forward", cdna[f:t])
                    optional_primers.append(tmp_primer)
                if th!=0.5:
                    if ((index+(l*(1-th))) >= 0) and ((index-(l*th)) < len(cdna)):
                        f = int(round((index - (l * (1 - th))), 0))
                        t = int(round((index + (l * th)), 0))
                        tmp_primer = Primer("forward", cdna[f:t])
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

    print len(primers)
    return primers

# ======================= tm test ===========================================


def tm_test(primer):
    if (primer.primer_tm() < 40) or (primer.primer_tm() > 70):
        return False
    return True

# ======================= %GC test ===========================================


def gc_test(primer):
    if (primer.precent_gc() < 40) or (primer.precent_gc() > 60):
        return False
    return True

# ======================= %syntaxTests test ===========================================


def syntax_tests(primer):
    if primer.sequnce[0] == 'G':
        return False

    if "GGGG" in primer.sequnce:
        return False

    if "CCCC" in primer.sequnce:
        return False

    if "AAAAAA" in primer.sequnce:
        return False

    consecutive_c = 0
    for index in range(len(primer.sequnce)):
        if primer.sequnce[index] == 'C':
            consecutive_c += 1
        else:
            consecutive_c = 0
        if consecutive_c == 2:
            return False    # primer cant have 2 consecutive C in the middle
    return True

# ======================= Primer score ============================================


def primers_score(primers):
    for primer in primers:
        # length score
        leng_avg = 20
        leng_stdev = 2.3
        leng_score_weight = 20
        if leng_avg-leng_stdev <= primer.length <= leng_avg+leng_stdev:
            primer.add_score(leng_score_weight)
        elif primer.length < leng_avg:
            primer.add_score(primer.length / leng_avg * leng_score_weight)
        else:
            primer.add_score((1-((primer.length / leng_avg )-1)) * leng_score_weight)

        # Tm score
        tm_avg = 57
        tm_stdev = 1.8
        tm_score_weight = 20
        if tm_avg-leng_stdev <= primer.primer_tm() <= tm_avg+tm_stdev:
            primer.add_score(tm_score_weight)
        elif primer.primer_tm() <= tm_avg:
            primer.add_score(primer.primer_tm() / tm_avg * tm_score_weight)
        else:
            primer.add_score((1 - ((primer.primer_tm() / tm_avg) - 1)) * tm_score_weight)


        # GC percent score
        gc_avg = 56
        gc_stdev = 8.3
        gc_score_weight = 20
        if gc_avg-leng_stdev <= primer.precent_gc() <= gc_avg+gc_stdev:
            primer.add_score(gc_score_weight)
        elif primer.precent_gc() <= gc_avg:
            primer.add_score(primer.precent_gc() / gc_avg * gc_score_weight)
        else:
            primer.add_score((1 - ((primer.precent_gc() / gc_avg) - 1)) * gc_score_weight)

# ======================= MAIN ====================================================


if __name__ == '__main__':
    #if len(sys.argv) == 3:
    #   species, symbol = sys.argv[1:]
    #else:
    #    species, symbol = 'human', 'BRAF'
    species = raw_input("Enter specie:")
    symbol = raw_input("Enter symbol:")
    transcript_data(species, symbol)