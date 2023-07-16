#!/usr/bin/env python3


from Bidding import *
from gameTreeVisualize import *
import sys


def run_test(north_secret_count, east_secret_count, chance_distribution, 
    par_table, max_par_value = None, visual = False):

    
    if max_par_value is None:
        max_par_value = north_secret_count
    else:
        max_par_value = min(max_par_value, north_secret_count)

    maxmin_det_value, maxmin_det_all_strategies = maxminDet(par_table, chance_distribution, north_secret_count, east_secret_count)
    
    maxmin_det_strategy_north = maxmin_det_all_strategies[0] # picking the 1st strategy in the list
    
    print("Maxmin value over deterministic strategies :", maxmin_det_value)
    print("----------------")

    print(maxmin_det_strategy_north)

    print("----------------")

    minmax_det_value, minamx_det_all_strategies = minmaxDet(par_table, chance_distribution, north_secret_count, east_secret_count)
    minmax_det_strategy_east = minamx_det_all_strategies[0] # picking the 1st strategy in the list

    
    print("Minmax value over deterministic strategies :", minmax_det_value)
    print("----------------")
    print(minmax_det_strategy_east)

    print("----------------")

    constraints = [('=',{0:{1:{'0':1}}},0)]
    game_value_beh, optimal_beh_strategy_north = behOpt(chance_distribution, par_table, north_secret_count, east_secret_count, 2, constraints)
    
    print("Game value : ", game_value_beh)
    print(optimal_beh_strategy_north)
    print("----------------")

    if visual:
        drawStrategyTree(north_secret_count, east_secret_count, chance_distribution, par_table, max_par_value,
         maxmin_det_strategy_north, minmax_det_strategy_east)


def parse_bidding_game_input(file_name): # returns a bidding game from a input file
    fh = open(file_name, 'r')
    north_secret_count = int(fh.readline().strip())
    east_secret_count =  int(fh.readline().strip())
    chance_dist = {}
    par_table = {}
    for line in fh:
        entry_list_of_string = line.strip().split()
        prob = float(entry_list_of_string.pop(2))

        [sec_N, sec_E, par_value_N, par_value_E] = map(int, entry_list_of_string)

        if not(0 <= sec_N < north_secret_count) or not(0 <= sec_E < east_secret_count):
            raise ValueError("Incorrect secret entry ", sec_N , sec_E)
        if (sec_N,sec_E) not in par_table:
            par_table[(sec_N,sec_E)] = [par_value_N, par_value_E]
        else:
            raise KeyError("Duplicate entry in par table", sec_N , sec_E)

        if 0 <= prob <= 1:
            chance_dist[(sec_N,sec_E)] = prob
        else:
            raise ValueError("Chance probability incorrect", prob)

    if sum(chance_dist.values()) != 1.0:
        raise ValueError("Chance probabilities do not add up to 1")
    if len(par_table) < north_secret_count*east_secret_count:
        raise KeyError("Missing values in par table")



    fh.close()
    max_par = 0
    for parlist in par_table.values():
        max_par = max(max_par, max(parlist))
    return (north_secret_count, east_secret_count,chance_dist, par_table, max_par)



if __name__ == "__main__":
    game = sys.argv[1]
    north_secret_count, east_secret_count,chance_dist, par_table, max_par = parse_bidding_game_input(game)
    run_test(north_secret_count, east_secret_count,chance_dist, par_table, max_par,visual=True)
