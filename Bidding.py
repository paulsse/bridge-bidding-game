#!/usr/bin/env python3


import itertools
from scipy.optimize import linprog
from fractions import Fraction


def StratSetN(north_secret_count, k): # generate set of all pure strategies of N ; k - highest val in par table
    Sub_strategy = {} # Sub_strategy[(P, b)] is set of all strategy subtrees of N,  where Player P will play after last bid b
    # P = 0 is N; P = 1 is E
    Sub_strategy[(0, k+1)] = [{'':0}]
    for i in range(k, -1, -1):
        statsE = []
        nextNstats = []
        for j in range(i+1, k+2): # E can play moves  i+1, .., k+1
            nextNstats.append(Sub_strategy[(0, j)]) #N-rooted subtrees after E plays j
        for statTuple in itertools.product(*nextNstats): #cross product of the subtrees
            branchE = {}
            for j in range( i+1, k+2):
                movej = statTuple[j- i-1]
                for hist in movej.keys():
                    branchE[str(j)+hist] = movej[hist]
            statsE.append(branchE)
        Sub_strategy[(1,  i)] = statsE
        statsN = [{'':k+1}]
        if i != 0:
            statsN.append({'':0})  
        # N passes or plays k+1
        for j in range( i+1, k+1): # N can play  i+1, ..., k
            for branchE in Sub_strategy[(1, j)]:
                temp = {'':j}
                for hist in branchE.keys():
                    temp[str(j)+hist] = branchE[hist]
                statsN.append(temp)
        if  i==0:
            for branchE in Sub_strategy[(1, 0)]:
                temp = {'':0}
                for hist in branchE.keys():
                    temp['0'+hist] = branchE[hist]
                statsN.append(temp)


        Sub_strategy[(0,  i)]= statsN
    final_strategy_Set = []
    combinations = itertools.product(*([Sub_strategy[(0, 0)]]*north_secret_count))
    for comb in combinations:
        st = {}
        for signal in range(north_secret_count):
            st[signal] = comb[signal]
        final_strategy_Set.append(st)
    return(final_strategy_Set)
        
def StratSetE(east_secret_count, k): # generate set of all pure strategies of E ; k - highest val in par table
    Sub_strategy = {} # Sub_strategy[(P, k)] is set of all strategy subtrees of E,  where Player P will play after last non-zero bid k
    # 0 is N; 1 is E
    Sub_strategy[(1, k+1)] = [{'':0}]
    for n in range(k, -1, -1):

        statsE = [{'':0}, {'':k+1}] # N passes or plays k+1
        for j in range(n+1, k+1): # N can play n+1, ..., k
            for branchN in Sub_strategy[(0, j)]:
                temp = {'':j}
                for hist in branchN.keys():
                    temp[str(j)+hist] = branchN[hist]
                statsE.append(temp)
        Sub_strategy[(1, n)]= statsE

        statsN = []
        nextEstats = []
        if n == 0:
            nextEstats.append(Sub_strategy[1, 0])
            for j in range(n+1, k+2): # N can play moves n+1, .., k+1
                nextEstats.append(Sub_strategy[(1, j)]) #E-rooted subtrees after N has played the move
            for statTuple in itertools.product(*nextEstats): #cross product of the subtrees
                branchE = {}
                for j in range(0, k+2):
                    movej = statTuple[j-n]
                    for hist in movej.keys():
                        branchE[str(j)+hist] = movej[hist]
                statsN.append(branchE)
            Sub_strategy[(0, n)] = statsN

            continue
        for j in range(n+1, k+2): # N can play moves n+1, .., k+1
                nextEstats.append(Sub_strategy[(1, j)]) #E-rooted subtrees after N has played the move
        for statTuple in itertools.product(*nextEstats): #cross product of the subtrees
            branchE = {}
            for j in range(n+1, k+2):
                movej = statTuple[j-n-1]
                for hist in movej.keys():
                    branchE[str(j)+hist] = movej[hist]
            statsN.append(branchE)
        Sub_strategy[(0, n)] = statsN
    final_strategy_Set = []
    combinations = itertools.product(*([Sub_strategy[(0, 0)]]*east_secret_count))
    for comb in combinations:
        st = {}
        for signal in range(east_secret_count):
            st[signal] = comb[signal]
        final_strategy_Set.append(st)
    return(final_strategy_Set)

def parPayoff(secN, secE, parTable, declarer, bid):
    par = parTable[(secN, secE)][declarer]
    if bid <= par:
        return(bid*(-1)**declarer)
    else:
        return((par-bid)*(-1)**declarer)


def bestResponseLocal(player,  opponent_strategy,  maxPar, parTable, chance, sec, history = ''): 
 #best response of player against opponent_strategy after observing sec and history
    secIndex = {} #indexing secrets w.r.t player
    secIndex[player] = sec
    if len(history) > 1 and history[-1] == '0': # Calculation of terminal payoff
        declarer = len(history)%2
        bid = int(history[-2]) #declarer's bid
        payoff = 0 
        for opponentSec in opponent_strategy.keys(): 
            secIndex[(1+player)%2] = opponentSec
            secN, secE = secIndex[0], secIndex[1]
            if declarer == player:
                opponentHist = history[:-1]
                opponentMove = 0
            else:
                opponentHist = history[:-2]
                opponentMove = int(history[-2])
            try:
                if opponent_strategy[opponentSec][opponentHist]==opponentMove: #checking if opponentonent's history is consistent
                    payoff+= chance[(secN, secE)]*parPayoff(secN, secE, parTable, declarer, bid)               
            except KeyError:
                continue
        return((payoff, []))
    elif len(history)%2 != player:
        nextHists = set()
        for opponentSec in opponent_strategy.keys():
            try:
                nextHists.add(history+str(opponent_strategy[opponentSec][history]))
            except KeyError:
                continue
        localVal = 0
        local_strategy_Set = []
        stratList = []
        for h in nextHists:
            opponentBid = int(h[-1])
            (val, localOpt_strategy) = bestResponseLocal(player, opponent_strategy, maxPar, parTable, chance, sec, h)
            localVal += val
            for st in localOpt_strategy:
                endHistList = list(st.keys())
                for endHist in endHistList:
                    st[str(opponentBid)+endHist] = st.pop(endHist)
            if len(localOpt_strategy)>0:
                stratList.append(localOpt_strategy)
        stratTupleSet = itertools.product(*stratList)
        for stratTuple in stratTupleSet:
            product_strategy = {}
            for st in stratTuple:
                for hist in st.keys():
                    product_strategy[hist] = st[hist]
            if len(product_strategy) > 0:
                local_strategy_Set.append(product_strategy)

        return(localVal, local_strategy_Set)
    else:
        if history == '':
            nextBid = 1
        else:
            nextBid = int(history[-1])+1
        (passPayoff, local_strategy_Set) = bestResponseLocal(player, opponent_strategy, maxPar, parTable, chance, sec, history + '0')
        if len(local_strategy_Set) > 0:
            for st in local_strategy_Set:
               endHistList = list(st.keys())
               for endHist in endHistList:
                   st['0'+endHist] = st.pop(endHist)
                   st[''] = 0
        else:
            local_strategy_Set.append({'':0})
        
        val = passPayoff
        for bid in range(nextBid, maxPar + 2): # player can bid between nextBid and maxPar +1
            (bestLocalVal, localOpt_strategy) =  bestResponseLocal(player, opponent_strategy, maxPar, parTable, chance, sec, history + str(bid))
            if (bestLocalVal - val)*(-1)**player > 0: #checking if local opt is better than previous based on parity of player
                val = bestLocalVal
                for st in localOpt_strategy:
                    endHistList = list(st.keys())
                    for endHist in endHistList:
                        st[str(bid)+endHist] = st.pop(endHist)
                        st[''] = bid
                if len(localOpt_strategy) == 0:
                    localOpt_strategy.append({'':bid})
                local_strategy_Set = localOpt_strategy

            elif bestLocalVal == val:
                for st in localOpt_strategy:
                    endHistList = list(st.keys())
                    for endHist in endHistList:
                        st[str(bid)+endHist] = st.pop(endHist)
                        st[''] = bid
                if len(localOpt_strategy) == 0:
                    localOpt_strategy.append({'':bid})
                local_strategy_Set.extend(localOpt_strategy)
        return(val, local_strategy_Set)

def bestResponse(player, opponent_strategy, maxPar, parTable, chance, secrets): 
    sub_strategy_List = [] 
    optVal = 0
    for sec in range(secrets):
        (localVal, local_strategy_Set) = bestResponseLocal(player,  opponent_strategy,  maxPar,  parTable,  chance,  sec)
        optVal += localVal
        sub_strategy_List.append(local_strategy_Set)
    opt_strategy_Comb = itertools.product(* sub_strategy_List)
    opt_strategy_Set = []
    for stratTuple in opt_strategy_Comb:
        stg = {}
        for sec in range(secrets):
            stg[sec] = stratTuple[sec]
        opt_strategy_Set.append(stg)
    return(optVal, opt_strategy_Set)


def maxminDet(par_table, chance, north_secret_count, east_secret_count): # returns maxmin value over pure stgs and set of maxmin optimal stgs
    k = 0  # k is the max par value in par table
    for v in par_table.values():
        k = max(k, max(v))

    Nstrategies = StratSetN(north_secret_count, k)
    l = len(Nstrategies)
    #print("Total strategies of N :", l)
    maxminpayoff = -(k+1)  # initializing maxmin value variable with lowest possible payoff
    maxmin_strategy_list = []
    i = 1
    for Nstrategy in Nstrategies:
        #print("Evaluating strategy " + str(i) + " of " + str(l) + " strategies" )
        brVal, brStg = bestResponse(1, Nstrategy, k, par_table, chance, east_secret_count)
        if brVal > maxminpayoff:
            maxminpayoff = brVal
            maxmin_strategy_list = [Nstrategy]
        elif maxminpayoff == brVal:
            maxmin_strategy_list.append(Nstrategy)
        #print("Maxmin Payoff is ", brVal)
        i+=1
    return(maxminpayoff, maxmin_strategy_list)

def minmaxDet(par_table, chance, north_secret_count, east_secret_count): # returns minmax value pure stgs and set of minmax optimal stgs
    k = 0  # k is the max par value in par table
    for v in par_table.values():
        k = max(k, max(v))

    Estrategies = StratSetE(east_secret_count, k)
    l = len(Estrategies)
    #print("Total strategies of E :", l)
    Nstrategies = StratSetN(north_secret_count, k)
    minmax_det_payoff = k+1  # initializing minmax value variable with highest possible payoff
    minmax_strategy_list = []
    i = 1
    for Estrategy in Estrategies:
        #print("Evaluating strategy " + str(i) + " of " + str(l) + " strategies" )
        brVal, brStgs = bestResponse(0, Estrategy, k, par_table, chance, north_secret_count)
        if brVal < minmax_det_payoff:
            minmax_det_payoff = brVal
            minmax_strategy_list = [Estrategy]
        elif minmax_det_payoff == brVal :
            minmax_strategy_list.append(Estrategy)
        #print("Minmax Payoff is ", brVal)
        i+=1
    return(minmax_det_payoff, minmax_strategy_list)



# Setup for game value calculation
def genSeqList(k): # increasing sequence generator
    if k == 0:
        return(['0'])
    else:
        subSeq = genSeqList(k-1)
        kSeq = []
        for seq in subSeq:
            kSeq.append(seq+str(k))
        subSeq.extend(kSeq)
        subSeq.append(str(k))
        return(subSeq) 

def hIndex(maxPar): #indexing histories in list
    nodeHists = genSeqList(maxPar + 1)
    Hists = []
    for hist in nodeHists:
        Hists.append(hist+'0')
    Hists.extend(nodeHists)
    Hists.append('')
    nh = len(Hists)
    Hists.sort()
    IndexDict = {}
    for i in range(nh):
        IndexDict[Hists[i]]=i
    return(Hists, IndexDict)


def coeffVector(coeffDict, histIndex, north_secret_count, east_secret_count): 
    #creating coefficient vector out of coefficient dictionary for LP constraints
    n = len(histIndex)
    vect = [0]*(north_secret_count+east_secret_count)*n
    for player in coeffDict.keys():
        for sec in coeffDict[player].keys():
            for hist in coeffDict[player][sec]:
                cf = coeffDict[player][sec][hist]
                vect[n*player*north_secret_count + n*sec + histIndex[hist]] = cf
    # (realization variables, value variables)
    return(vect) 

def isTerminal(hist): #returns true if hist is terminal history
    return(len(hist) > 1 and int(hist[-1]) == 0)

def nextBids(hist, maxPar): # returns list of next playable bids after hist
    if isTerminal(hist):
        return([])
    else:
        if len(hist) == 0:
            bids = list(range(maxPar+2))
        else:
            lastBid = int(hist[-1])
            bids = []
            bids.append(0)
            for i in range(lastBid+1, maxPar+2):
                bids.append(i)
        return(bids) # list of next possible bids

def behOpt(chance, parTable, north_secret_count, east_secret_count, maxPar, startVals = []): #
    # returns value of the game and optimal behavioral strategy
    Hists, histIndex = hIndex(maxPar)
    IneqVect = []
    IneqConst = []
    EqVect = []
    EqConst = []
    for hist in Hists:
        if isTerminal(hist):
            #terminal equations
            for secE in range(east_secret_count):
                Cf = {1:{secE:{hist:-1}}}
                Cf0 = {}
                for secN in range(north_secret_count):
                    bid = int(hist[-2])
                    declarer = len(hist)%2
                    normalizedTerminalPayoff = chance[(secN, secE)]*parPayoff(secN, secE, parTable, declarer, bid)
                    hisCf = {hist:chance[(secN, secE)]*parPayoff(secN, secE, parTable, declarer, bid)}
                    Cf0[secN] = hisCf
                Cf[0] = Cf0
                EqVect.append(coeffVector(Cf, histIndex, north_secret_count, east_secret_count))
                EqConst.append(0)
        else:
            #realization plans
            for sec in range(north_secret_count):
                if len(hist)==0:               
                    cfEq = {0:{sec:{'':1}}}
                    EqVect.append(coeffVector(cfEq, histIndex, north_secret_count, east_secret_count))
                    EqConst.append(1)

                if len(hist)%2 == 0: 
                    #sum equality to 1
                    bids = nextBids(hist, maxPar)
                    hisCf = {}
                    hisCf[hist] = 1
                    for b in bids:
                        hisCf[hist+str(b)] = -1
                    Cf = {0:{sec:hisCf}}
                    EqVect.append( coeffVector(Cf, histIndex, north_secret_count, east_secret_count))
                    EqConst.append(0)
                else:
                    #plan propagation
                    bids = nextBids(hist, maxPar)
                    for b in bids:
                        hisCf = {}
                        hisCf[hist] = 1
                        hisCf[hist+str(b)] = -1
                        Cf = {0:{sec:hisCf}}
                        EqVect.append( coeffVector(Cf, histIndex, north_secret_count, east_secret_count))
                        EqConst.append(0)
            #valuation equalities
            if len(hist)%2 == 0:
                for sec in range(east_secret_count):
                    bids = nextBids(hist, maxPar)
                    hisCf = {}
                    hisCf[hist] = 1
                    for b in bids:
                        hisCf[hist+str(b)] = -1
                    Cf = {1:{sec:hisCf}}
                    EqVect.append( coeffVector(Cf, histIndex, north_secret_count, east_secret_count))
                    EqConst.append(0)

            else:
                # value local optimality inequalities
                for sec in range(east_secret_count):
                    bids = nextBids(hist, maxPar)
                    for b in bids:
                        hisCf = {}  
                        hisCf[hist] = 1
                        hisCf[hist+str(b)] = -1
                        Cf = {1:{sec:hisCf}}
                        IneqVect.append( coeffVector(Cf, histIndex, north_secret_count, east_secret_count))
                        IneqConst.append(0)

    Cf = {}
    objCf = {}

    #Additional constraints for initial values
    for (sign, Cf, const) in startVals:
        if sign == '=':
            EqVect.append(coeffVector(Cf, histIndex, north_secret_count, east_secret_count))
            EqConst.append(const)
        else:
            IneqVect.append(coeffVector(Cf, histIndex, north_secret_count, east_secret_count))
            IneqConst.append(const)
    for secE in range(east_secret_count):
        objCf[secE] = {'':-1}
    Cf[1] = objCf
    # objective vector
    Obj = coeffVector(Cf, histIndex, north_secret_count, east_secret_count)
    n = len(histIndex)

    # variable bounds
    Bounds =  [(0, 1)]*(n*north_secret_count) 
    vals = [(None, None)]*(n*east_secret_count)
    Bounds.extend(vals)

    gameOpt = linprog(Obj, A_ub = IneqVect,  b_ub = IneqConst, 
        A_eq=EqVect, b_eq= EqConst, bounds = Bounds, method = "highs")

    optStg = {}
    soln = gameOpt.x
    
    n = len(Hists)
    for secN in range(north_secret_count):
        optStg[secN] = {}
        for h in Hists:
            totalProb = soln[n*secN + histIndex[h]]
            if len(h)%2 == 0 and not isTerminal(h) and (totalProb > 0):
                optStg[secN][h] = {}
                for b in nextBids(h, maxPar):
                    prob = soln[n*secN + histIndex[h+str(b)]]/totalProb
                    optStg[secN][h][b] = prob
    optVal = -gameOpt.fun
    return(optVal, optStg) # returns game value profile




