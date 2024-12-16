'''
Two State HMM Implementations
(n) Seth Temple (e) sdtemple at uw dot edu
'''

from math import log2
from copy import copy    

def log_emissions(e):
    '''
    Convert emission probs to log scale
    :param e: emission probs
    :return: log emission probs
    :rtype: dict
    '''
    return {key: log2(val) for (key, val) in e.items()}

def log_transitions(t):
    '''
    Convert transition probs to log scale
    :param t: transition probs
    :return: log transition probs
    :rtype: list
    '''
    T = []
    for tt in t:
        T.append((log2(tt[0]), log2(tt[1])))
    return T

# not optimal runtime; take a bathroom break
def viterbi_parse(obs, i, t, e):
    '''
    Determine the Viterbi parse for a sequence
    :param obs: observation sequence
    :param i: initiation probs
    :param t: transition probs
    :return:
    :rtype:
    '''
    top = [0]
    tscore = i[0] + e[0][obs[0][1]]
    bot = [1]
    bscore = i[1] + e[1][obs[0][1]]
    for i in range(1, len(obs)):
        candidate00 = tscore + t[0][0] + e[0][obs[i][1]]
        candidate01 = tscore + t[0][1] + e[1][obs[i][1]]
        candidate10 = bscore + t[1][0] + e[0][obs[i][1]]
        candidate11 = bscore + t[1][1] + e[1][obs[i][1]]

        ctop = copy(top)
        cbot = copy(bot)
        if candidate00 >= candidate10:
            tscore = candidate00
            top.append(0)
        else:
            tscore = candidate10
            top = cbot
            top.append(0)

        if candidate01 >= candidate11:
            bscore = candidate01
            bot = ctop
            bot.append(1)
        else:
            bscore = candidate11
            bot.append(1)

        '''
        if (i % 10000) == 0:
            print(i)
        '''

    if tscore >= bscore:
        return (tscore, top)
    else:
        return (bscore, bot)

def hist_states(parse):
    '''
    Create a histogram for # times in each state
    :param parse: state sequence of Viterbi parse
    :return: histogram
    :rtype: dict
    '''
    d = {0:0, 1:0}
    for i in range(len(parse)):
        d[parse[i]] += 1
    return d

def hist_segments(parse, obs):
    '''
    Create a histogram for # state segments
    :param parse: state sequence of Viterbi parse
    :param obs: observation sequence
    :return: histogram and state 1 segments
    :rtype: tuple(dict, list)
    '''
    d = {0:0, 1:0}
    s = []
    pstate = parse[0]
    ppos = obs[0][0]
    for i in range(1, len(parse)):
        cstate = parse[i]
        if ((pstate == 0) and (cstate == 1)):
            d[0] += 1
            ppos = obs[i][0]
            pstate = cstate
        elif ((pstate == 1) and (cstate == 0)):
            d[1] += 1
            cpos = obs[i - 1][0]
            s.append((ppos, cpos))
            pstate = cstate

    # sort segments by length
    def seg_sort(tup):
        return tup[1] - tup[0]
    s.sort(key = seg_sort, reverse = True)
    
    return (d, s)
