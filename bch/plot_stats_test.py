import sys
import random

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from bch import *
from simulate import *

import os
import functools

def get_random_messages(count, n):
    messages = []
    for message in range(count):
        candidate = random.getrandbits(n)

        # is_present = False
        # for test in messages:
        #     if candidate == test:
        #         is_present = True
        #         break
        
        # if is_present:
        #     message -= 1
        #     continue
        # else:
        messages.append(candidate)
    
    return messages

n = 15
for t in range(2, 10):
    bch = BCH(n, 2 * t + 1, 1, get_primitive_polynomial(4))
    
    # encoding
    blocks = []
    encoded_codewords = []
    for block in get_random_messages(count=100, n=bch.k):
        blocks.append(block)
        codeword = bch.encode(block)
        if VERBOSE:
            print("Block: {0:>0{a}b}, codeword: {1:>0{b}b}".format(block, codeword, a=bch.n - bch.t * bch.power, b=bch.n))
        encoded_codewords.append(codeword)

    # plot
    mpl.rc('axes', axisbelow=True)
    fig, ax = plt.subplots(1, 1, figsize=(9.6, 7.5))
    for ind, spine in ax.spines.items():
        if ind != 'top':
            spine.set_visible(False)
    ax.set_adjustable('box')
    ax.invert_yaxis()

    ax.set_xlim(0, len(encoded_codewords))
    ax.set_ylim(1 - .5, n + .5)

    ax.xaxis.set_major_locator(ticker.LinearLocator(10))
    ax.xaxis.set_major_formatter(ticker.PercentFormatter(xmax=len(encoded_codewords), symbol=''))
    ax.xaxis.set_label_position('top')

    ax.set_yticks([i for i in range(1, n + 1)])
    ax.tick_params(top=True, bottom=False, left=False, right=False,
                labeltop=True, labelbottom=False, labelleft=True)
    ax.tick_params(which='minor', bottom=False, top=True)

    # ax.legend(ncol=3, loc='upper center', bbox_to_anchor=(.5, 0))

    ax.set_xlabel('Decoding Results, %')
    ax.set_ylabel('Errors Made')

    for error_number in range(t, t + 2):
        print("Running with t = {}".format(error_number))
        # distortion
        distorted_codewords = []
        for codeword in encoded_codewords:
            error_vector = get_random_number_of_hamming_weight(length=bch.n, weight=error_number)
            distorted_codeword = codeword ^ error_vector

            if VERBOSE:
                print("Codeword: {0:0>{a}b}, error vector: {1:0>{a}b}, distorted codeword: {2:0>{a}b}".
                    format(codeword, error_vector, distorted_codeword, a=bch.n))
            distorted_codewords.append(distorted_codeword)
        
        # correcting
        right  = 0
        wrong  = 0
        cancel = 0
        for codeword in range(len(distorted_codewords)):
            try:
                message = bch.decode(distorted_codewords[codeword])
            except Exception as e:
                cancel += 1
                continue

            if VERBOSE:
                print("Codeword: {0:0>{a}b}, Decoded message: {1:0>{a}b}".
                    format(blocks[codeword], message, a=bch.k))
            if blocks[codeword] == message:
                right += 1
            else:
                wrong += 1
                message = bch.decode(distorted_codewords[codeword])

        assert right + wrong + cancel == len(encoded_codewords)

        ax.barh(error_number, left=0, width=right,
                height=.9, color='C2', label='Right')
        ax.barh(error_number, left=right, width=wrong,
                height=.9, color='C3', label='Wrong')
        ax.barh(error_number, left=right + wrong, width=cancel,
                height=.9, color='C0', label='Cancel')

    fig.savefig('images/statistics_bch-n-{}-k-{}-t-{}-prim-{}.png'.format(bch.n, bch.k, bch.t, bch.power),
                bbox_inches='tight', dpi=200)




