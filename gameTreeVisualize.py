#!/usr/bin/env python3


from Bidding import *
from graphviz import Digraph

alph = {0:{0:'a',1:'b',2:'c',3:'d'},1:{0:'x',1:'y',2:'z',3:'w'}}


# generates the game tree and color the edges conforming to the strategies provided

def drawStrategyTree(north_secret_count, east_secret_count, chance_distribution, par_table,
                         max_par_value, strategy_North, strategy_East,  name = 'GameStrategyTree'):
    game = Digraph()
    game.node("chance",shape = "triangle") # Chance node
    for (secN,secE) in chance_distribution.keys():
        game.edge("chance",alph[0][secN]+ " "+alph[1][secE],label = str(chance_distribution[(secN,secE)]))

    all_histories = hIndex(max_par_value)[0]
    for (secN,secE) in chance_distribution.keys():
        secStr = alph[0][secN]+ " "+alph[1][secE]
        for hist in all_histories:

            if isTerminal(hist): #leaf nodes
                bid = int(hist[-2])
                declarer = len(hist)%2
                game.node(secStr + hist,shape = "point", fontcolor = "darkgreen", xlabel = str(parPayoff(secN,secE,par_table,declarer,bid)))

            elif len(hist)%2 == 0: # North node
                game.node(secStr + hist,label = secStr +'\n'+ hist,shape = "circle",style='filled',color='lightblue')
            else: # East node
                game.node(secStr + hist,label = secStr + '\n'+ hist,shape = "square",style='filled', color='lightgrey')
            
            # Coloring Edges
            if len(hist) > 0:
                edge_color = "black"
                font_color = "black"
                if len(hist)%2 == 0: #last bid was played by East
                    if hist[:-1] in strategy_East[secE] and strategy_East[secE][hist[:-1]] == int(hist[-1]): # last bid was played using stratgy of East
                        edge_color = "red"
                        font_color = "red"

                if len(hist)%2 == 1: #last bid was played by North
                    if hist[:-1] in strategy_North[secN] and strategy_North[secN][hist[:-1]] == int(hist[-1]): # last bid was played using stratgy of North
                        edge_color = "red"
                        font_color = "red"


                game.edge(secStr + hist[:-1],secStr + hist,label = hist[-1],color = edge_color, fontcolor = font_color)
    
    
    game.render(name, view = True)  
