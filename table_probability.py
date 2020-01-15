#!/usr/bin/env python3 
#-*- coding: utf-8 -*-

import sys
import random

from bch import *

import os
import functools

def get_random_messages(count, n):
    messages = []
    for _ in range(count):
        candidate = random.getrandbits(n)
        messages.append(candidate)
    
    return messages

def get_hamming_weight(num):
    length = len(bin(num)) - 2
    weight = 0

    for i in range(length):
        if num & 1 == 1:
            weight += 1
        num >>= 1

    return weight

def error_probability():
    print("│     m     │    bch    │     t     │    t + 1   │ ")
    print("├───────────┼───────────┼───────────┼────────────┤ ")

    procent = []
    for m in range(4, 11):

        n = 31
        t = 0
        if (m == 4):
            n = 15
            t = 2
        elif (m >= 5 and m <= 6):
            t = 3
        elif (m >= 7):
            t = 2

        bch = BCH(n, 2 * t + 1, 1, get_primitive_polynomial(m))

        # encoding
        blocks = []
        encoded_codewords = []
        for block in get_random_messages(count=10000, n=bch.k):
            blocks.append(block)
            codeword = bch.encode(block)
            if VERBOSE:
                print("Block: {0:>0{a}b}, codeword: {1:>0{b}b}".format(block, codeword, a=bch.n - bch.t * bch.power, b=bch.n))
            encoded_codewords.append(codeword)

        for error_number in range(t, t + 2):
            if VERBOSE:
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

            bits_corrupted = 0
            bits_all = 0

            for codeword in range(len(distorted_codewords)):
                try:
                    message, decoded_codeword = bch.decode_full(distorted_codewords[codeword])
                except Exception as e:
                    cancel += 1
                    continue

                if VERBOSE:
                    print("Codeword: {0:0>{a}b}, Decoded message: {1:0>{a}b}".
                        format(blocks[codeword], message, a=bch.k))
                
                bits_all += bch.n
                if blocks[codeword] == message:
                    right += 1
                else:
                    wrong += 1
                    bits_corrupted += get_hamming_weight(distorted_codewords[codeword] ^ decoded_codeword)

            assert right + wrong + cancel == len(encoded_codewords)

            procent.append(wrong)

        print("│     {:<2}    │  {:<2} {:<2} {:<2} │     {:<3}   │    {:<4}    │ ".format(
            bch.power,
            bch.n, bch.k, bch.dist,
            procent[0],
            procent[1])
        )

        procent = []

if __name__ == '__main__':
    error_probability()
