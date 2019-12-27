#!/usr/bin/env python3 
#-*- coding: utf-8 -*-
import math
import random

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker

from bch import *

def get_random_messages(count, n):
    messages = []
    for _ in range(count):
        candidate = random.getrandbits(n)
        messages.append(candidate)
    
    return messages

def bits_to_signal(bits):
    signal_array = []
    for i in range(msb(bits), -1, -1):
        signal_array.append(1. if is_bit_set(bits, i) else -1.)
    return signal_array

def add_noise_to_signal(signal, sigma):
    return [random.gauss(i, sigma) for i in signal]

def signal_to_bits(signal):
    res = 0
    for i in signal:
        res <<= 1
        res ^= 1 if i > 0 else 0
    return res

def awgn(n_count, E_b_N_0):
    bch = BCH(512, 63, 1, get_primitive_polynomial(9))
    # bch = BCH(15, 5, 1, get_primitive_polynomial(4))
    
    code_speed = bch.k / (bch.n * 1.) 
    Eb = 1 / code_speed
    N_0 = Eb / E_b_N_0 # == 2 * sigma * sigma
    sigma = math.sqrt(N_0 / 2)

    counter = 0
    # VERBOSE = True
    for block in get_random_messages(count=n_count, n=bch.k):
        if VERBOSE:
            print("-----------------------------")
            print(bin(block)[2:])

        # encode
        codeword = bch.encode(block)
        if VERBOSE:
            print(bin(codeword)[2:])

        # create signal
        signal = bits_to_signal(codeword)
        if VERBOSE:
            print(signal)

        # noise
        distorted_signal = add_noise_to_signal(signal, sigma)
        if VERBOSE:
            print(distorted_signal)

        # recieve
        recv = signal_to_bits(distorted_signal)
        if VERBOSE:
            print(bin(recv)[2:])

        # decode
        message = bch.decode(recv)
        if VERBOSE:
            print(bin(message)[2:])

        if block != message:
            print("Error on iteration {}".format(counter))
        counter += 1

if __name__ == '__main__':
    awgn(n_count=5000, E_b_N_0=4)
