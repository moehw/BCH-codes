#!/usr/bin/env python3 
#-*- coding: utf-8 -*-

import math
import random
import multiprocessing

PRE_DEFINED = False
PLOT = False

if PLOT:
    import matplotlib.pyplot as plt

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

def get_sigma_by_SNR_and_speed(signal_to_noise_ratio, code_speed):
    return math.sqrt(1.0 / (2 * code_speed * math.pow(10, x / 10.0)))

def awgn_with_coding(iteration_count, signal_to_noise_ratio, n, bch=None):
    if bch is not None:
        code_speed = bch.k / (bch.n * 1.) 
        k = bch.k
    else:
        k = n
        code_speed = 1

    sigma = get_sigma_by_SNR_and_speed(signal_to_noise_ratio, code_speed)

    counter = 0
    errors_message = 0
    errors_word = 0

    for block in get_random_messages(count=iteration_count, n=k):
        if VERBOSE:
            print("-----------------------------")
            print(bin(block)[2:])

        # encode
        if bch is not None:
            codeword = bch.encode(block)
            if VERBOSE:
                print(bin(codeword)[2:])
        else:
            codeword = block

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
        if bch is not None:
            message, decoded_codeword = bch.decode_ex(recv)
            if VERBOSE:
                print(bin(message)[2:])

            if codeword != decoded_codeword:
                if VERBOSE:
                    print("Error word on iteration {}".format(counter))
                errors_word += 1

            if block != message:
                if VERBOSE:
                    print("Error message on iteration {}".format(counter))
                errors_message += 1
        else:
            if recv != block:
                if VERBOSE:
                    print("Error word on iteration {}".format(counter))
                errors_word += 1

        print("SNR: {}\t-> iteration: {}".format(signal_to_noise_ratio, counter))
        counter += 1

    return errors_word, errors_message

def awgn_worker(signal_to_noise_ratio, test_number, n, bch, return_dict):
    no_coding_result = (awgn_with_coding(iteration_count=test_number, signal_to_noise_ratio=signal_to_noise_ratio, n=n)[0] / (test_number * 1.))

    signal_to_noise_ratio_result = (awgn_with_coding(iteration_count=test_number, signal_to_noise_ratio=signal_to_noise_ratio, n=n, bch=bch)[0] / (test_number * 1.))

    return_dict[signal_to_noise_ratio] = no_coding_result, signal_to_noise_ratio_result

if __name__ == '__main__':
    x_array = [_ for _ in range(1, 14)]

    signal_to_noise_ratio_result = []
    no_coding_result = []

    if not PRE_DEFINED:
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        jobs = []

        test_number = 10

        n = 512
        bch = BCH(n, 63, 1, get_primitive_polynomial(9))
        
        for x in x_array:
            p = multiprocessing.Process(target=awgn_worker, args=(x, test_number, n, bch, return_dict))
            jobs.append(p)
            p.start()

        for proc in jobs:
            proc.join()
        
        return_list = sorted(return_dict.items())

        for key, value in return_list:
            no_coding_result.append(value[0])
            signal_to_noise_ratio_result.append(value[1])

        print("No coding: {}".format(no_coding_result))
        print("SNR: {}".format(signal_to_noise_ratio_result))

    else: #                             1    2    3     4     5     6     7     8     9     10    11    12   13
        signal_to_noise_ratio_result = [1.0, 1.0, 0.99, 0.42, 0.00, 0.01, 0.01, 0.0 , 0.0 , 0.00, 0.0 , 0.0, 0.0]
        no_coding_result             = [1.0, 1.0, 1.0 , 1.0 , 1.0 , 0.99, 0.97, 0.78, 0.44, 0.19, 0.05, 0.0, 0.0]

    if PLOT:
        plt.plot(x_array[0:len(signal_to_noise_ratio_result)], signal_to_noise_ratio_result, 'o') # points
        plt.plot(x_array[0:len(no_coding_result)]            , no_coding_result, '+')             # points
        plt.plot(x_array[0:len(signal_to_noise_ratio_result)], signal_to_noise_ratio_result, '-') # line
        plt.plot(x_array[0:len(no_coding_result)]            , no_coding_result, '--')            # line

        plt.ylabel('Probability of error')
        plt.xlabel('Eb/N0, dB')

        plt.yscale('log')
        plt.grid(True)

        plt.xticks(x_array[1::2])

        plt.show()
