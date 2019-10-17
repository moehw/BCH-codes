#!/usr/bin/env python2
from galois import *
from coding import *

i = {
    '0': 0,
    '1': 1,
    '2': 2,
    '3': 3,
    '4': 4,
    '5': 5,
    '6': 6,
    '7': 7,
    '8': 8,
    '9': 9,
    'a': 10,
    'b': 11,
    'c': 12,
    'd': 13,
    'e': 14,
    'f': 15,
    'g': 16,
    'h': 17,
    'i': 18,
    'j': 19,
    'k': 20,
    'l': 21,
    'm': 22,
    'n': 23,
}

to_letter = {}
for key, value in i.items():
    to_letter[value] = key


def print_letter(s):
    res = ""
    for num in s:
        try:
            res += to_letter[num]
        except KeyError:
            res += str(num)
        res += " "
    return res

def goley(g, y, replace_table, IS):

    # IS - information set - usually count of rows
    g = Matrix(data=g).to_GF(2)
    y = Matrix(data=y).to_GF(2)

    print("G:")
    print(g)
    print

    print("y:")
    print(y)
    print

    dist_list = []

    best_canonical = []
    best_replace   = []
    best_u         = []
    best_y_r       = []
    best_y2        = []
    best_dist      = 23

    # for every replace string
    for r in replace_table:

        print("--------------------------------------------------")

        g_replaced = g.copy()
        y_replaced = y.copy()

        for el in range(g_replaced.cols):
            g_replaced.copy_col(el, g, r[el])
            y_replaced.copy_col(el, y, r[el])
        
        # now we have replaced g and y with one replace-string

        print("G repl:")
        print(g_replaced)
        print

        # canonical form
        g_canonical = g_replaced.get_reduced_echelon()
        print("G canonical:")
        print(g_canonical)
        print

        u_replaced = y_replaced.submatrix(0, 0, 1, IS) # information set
        y2 = u_replaced * g_canonical
        print("Repl. string: " + print_letter(r))
        print("u  replaced : {0}".format(u_replaced))
        print("y' replaced : {0}".format(y_replaced))
        print("y'' (u*G)   : {0}".format(y2))

        dist = y2.distance(y_replaced)
        print("Distance H  : len {0}".format(len(dist)))
        dist_str = ""
        for x in dist:
            dist_str += to_letter.get(int(x[1])) + " "
        print dist_str
        print

        if best_dist > len(dist):
            best_canonical = g_canonical.copy()
            best_replace   = r
            best_u         = u_replaced.copy()
            best_y_r       = y_replaced.copy()
            best_y2        = y2.copy()
            best_dist      = len(dist)

        dist_list.append(len(dist))

    print("##################################################")
    print("All distances: {0}".format(dist_list))
    print
    print("The best one ({0}):".format(best_dist))
    print("G canonical:")
    print(best_canonical)
    print
    print("Repl. string: " + print_letter(best_replace))
    print("u  replaced : {0}".format(best_u))
    print("y' replaced : {0}".format(best_y_r))
    print("y'' (u*G)   : {0}".format(best_y2))
    dist = best_y2.distance(best_y_r)
    print("Distance H  : len {0} : {1}".format(len(dist), dist))

def test():
    
    g = [
        [1,0,0,0,0,1,1],
        [0,1,0,0,1,0,1],
        [0,0,1,0,1,1,0],
        [0,0,0,1,1,1,1]
    ]

    y = [[0,1,1,0,1,0,0]]

    repl = [[3-1,5-1,6-1,7-1,1-1,2-1,4-1]]

    goley(g, y, repl, 4)

def task1():
    # 129 variant

    #    0  1  2  3  4  5  6  7  8  9  a  b  c  d  e  f  g  h  i  j  k  l  m  n 
    y = [
        [0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0]
    ]

    #    0  1  2  3  4  5  6  7  8  9  a  b  c  d  e  f  g  h  i  j  k  l  m  n 
    g = [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1],
    ]

    replace_table = [
        [i['0'], i['1'], i['2'], i['3'], i['4'], i['5'], i['6'], i['7'], i['8'], i['9'], i['a'], i['b'], i['c'], i['d'], i['e'], i['f'], i['g'], i['h'], i['i'], i['j'], i['k'], i['l'], i['m'], i['n']],
        [i['0'], i['1'], i['2'], i['c'], i['d'], i['e'], i['3'], i['4'], i['5'], i['f'], i['g'], i['h'], i['i'], i['j'], i['k'], i['6'], i['7'], i['8'], i['l'], i['m'], i['n'], i['9'], i['a'], i['b']],
        [i['0'], i['1'], i['2'], i['i'], i['j'], i['k'], i['c'], i['d'], i['e'], i['6'], i['7'], i['8'], i['l'], i['m'], i['n'], i['3'], i['4'], i['5'], i['9'], i['a'], i['b'], i['f'], i['g'], i['h']],
        [i['0'], i['1'], i['2'], i['l'], i['m'], i['n'], i['i'], i['j'], i['k'], i['3'], i['4'], i['5'], i['9'], i['a'], i['b'], i['c'], i['d'], i['e'], i['f'], i['g'], i['h'], i['6'], i['7'], i['8']],
        [i['0'], i['1'], i['2'], i['9'], i['a'], i['b'], i['l'], i['m'], i['n'], i['c'], i['d'], i['e'], i['f'], i['g'], i['h'], i['i'], i['j'], i['k'], i['6'], i['7'], i['8'], i['3'], i['4'], i['5']],
        [i['0'], i['1'], i['2'], i['f'], i['g'], i['h'], i['9'], i['a'], i['b'], i['i'], i['j'], i['k'], i['6'], i['7'], i['8'], i['l'], i['m'], i['n'], i['3'], i['4'], i['5'], i['c'], i['d'], i['e']],
        [i['0'], i['1'], i['2'], i['6'], i['7'], i['8'], i['f'], i['g'], i['h'], i['l'], i['m'], i['n'], i['3'], i['4'], i['5'], i['9'], i['a'], i['b'], i['c'], i['d'], i['e'], i['i'], i['j'], i['k']],
        [i['c'], i['d'], i['e'], i['f'], i['g'], i['h'], i['i'], i['j'], i['k'], i['l'], i['m'], i['n'], i['0'], i['1'], i['2'], i['3'], i['4'], i['5'], i['6'], i['7'], i['8'], i['9'], i['a'], i['b']],
        [i['i'], i['j'], i['k'], i['6'], i['7'], i['8'], i['l'], i['m'], i['n'], i['9'], i['a'], i['b'], i['0'], i['1'], i['2'], i['c'], i['d'], i['e'], i['3'], i['4'], i['5'], i['f'], i['g'], i['h']],
        [i['l'], i['m'], i['n'], i['3'], i['4'], i['5'], i['9'], i['a'], i['b'], i['f'], i['g'], i['h'], i['0'], i['1'], i['2'], i['i'], i['j'], i['k'], i['c'], i['d'], i['e'], i['6'], i['7'], i['8']],
        [i['9'], i['a'], i['b'], i['c'], i['d'], i['e'], i['f'], i['g'], i['h'], i['6'], i['7'], i['8'], i['0'], i['1'], i['2'], i['l'], i['m'], i['n'], i['i'], i['j'], i['k'], i['3'], i['4'], i['5']],
        [i['f'], i['g'], i['h'], i['i'], i['j'], i['k'], i['6'], i['7'], i['8'], i['3'], i['4'], i['5'], i['0'], i['1'], i['2'], i['9'], i['a'], i['b'], i['l'], i['m'], i['n'], i['c'], i['d'], i['e']],
        [i['6'], i['7'], i['8'], i['l'], i['m'], i['n'], i['3'], i['4'], i['5'], i['c'], i['d'], i['e'], i['0'], i['1'], i['2'], i['f'], i['g'], i['h'], i['9'], i['a'], i['b'], i['i'], i['j'], i['k']],
        [i['3'], i['4'], i['5'], i['9'], i['a'], i['b'], i['c'], i['d'], i['e'], i['i'], i['j'], i['k'], i['0'], i['1'], i['2'], i['6'], i['7'], i['8'], i['f'], i['g'], i['h'], i['l'], i['m'], i['n']],
    ]

    goley(g, y, replace_table, 12)


if __name__ == '__main__':
    task1()
    #test()

