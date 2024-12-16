'''
.gbff Parsing for Site Models
(n) Seth Temple (e) sdtemple at uw dot edu
'''

from gbff_parse import *
from gsci_utility import *
from math import floor

# input
fn = '../D/Gallus_gallus_chr26.gbff.txt'
with open(fn, 'r') as f:
    to_features(f)
    to_origin(f)
    gnm = get_gnm(f)

# do math
cgnm = complement_gnm(gnm)
bctr = nucl_cts(gnm)
bfrq = nucl_freqs(bctr)
tctr = count_sites(fstrand, cstrand, gnm, cgnm)
tfrq = [freq_site(site) for site in tctr]
wts = [wt_site(site, bfrq) for site in tfrq]
max_scores = [max(d.values()) for d in wts]
cds_scores = [score_site(site, gnm, wts) for site in fstrand] + [score_site(site, cgnm, wts) for site in cstrand]
rcds_scores = [floor(scr) for scr in cds_scores]
cds_hist = {}
for scr in rcds_scores:
    try:
        cds_hist[scr] += 1
    except KeyError:
        cds_hist[scr] = 1

fpos = [site[10] for site in fstrand]
cpos = [site[10] for site in cstrand]

all_hist = {}
pos_list = []
lt = 0
n = len(gnm)
for i in range(10, n - 10):
    # forward
    site = [i + k for k in range(-10,11)]
    scr = score_site(site, gnm, wts)
    if scr >= 10:
        if site[10] not in fpos:
            pos_list.append((site[10] + 1, 0, scr))

    fscr = floor(scr)
    try:
        if fscr < -50:
            lt += 1
        else:
            all_hist[fscr] += 1
    except KeyError:
        all_hist[fscr] = 1

    # complement
    csite = [i + k for k in range(10,-11,-1)]
    cscr = score_site(csite, cgnm, wts)
    if cscr >= 10:
        if csite[10] not in cpos:
            pos_list.append((csite[10] + 1, 1, cscr))
    
    fcscr = floor(cscr)
    try:
        if fcscr < -50:
            lt += 1
        else:
            all_hist[fcscr] += 1
    except KeyError:
        all_hist[fcscr] = 1






