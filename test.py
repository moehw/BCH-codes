#!/usr/bin/env python3 
#-*- coding: utf-8 -*-

import sys

from bch import *

# Cases

def lab():
        primpoly = 0b1101101
        n = 63
        d = 11
        t = (d - 1) >> 1
        b = 55
        bch = BCH(n, d, b, primpoly)
        return bch

def lecture():
    bch = BCH(15, 5, 1, 0b10011)
    return bch

def book():
    generator = 0b10100110111
    n = 15
    d = 7
    bch = BCH(n, d, 1, 0b10011)
    return bch

# Tests

def test_bch():
    
    bchs     = [lab()        , book(), lecture()]
    messages = [0b11010000101, 0     , 0b1010101]

    for sample in range(len(bchs)):

        print("\nTest BCH №{}".format(sample + 1))
        bch = bchs[sample]
        message = messages[sample]

        print("Input: {} ({:b})".format(binary_to_string(message), message))

        result = bch.encode(message)
        print("Result encode: {} ({:b})".format(binary_to_string(result), result))

        # result = 0b10100100
        result ^= 0b10001001

        print("Distorted code: {} ({:b})".format(binary_to_string(result), result))

        result = bch.decode(result)
        # result = bch.bm(result)
        print("Result decode: {} ({:b})".format(binary_to_string(result), result))

        assert message == result, \
            "Message was encoded or decoded incorrect"

def test_berlekamp_massey_decode():
        print("\nTest Berlekamp-Massey")
        
        power = 4
        t = 2
        primitive_polynomial = 0b10011
        received_message = 0b000011001100011
        cyclotomic_cosets = get_cyclotomic_cosets(power=power)
        
        logarithm_table = get_logarithm_table(power=power, primitive_polynomial=primitive_polynomial)
        syndromes, is_error = get_syndromes(
            primitive_polynomial=primitive_polynomial,
            received_message=received_message,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithm_table=logarithm_table,
            power=power,
            t=t,
            b=1
        )
        if is_error:
            result = berlekamp_massey_decode(
                syndromes=syndromes,
                logarithm_table=logarithm_table,
                power=power,
                t=t,
                primitive_polynomial=primitive_polynomial)

        assert result == [0, 2, 14, -1]

# Testing

if __name__ == '__main__':
    test_berlekamp_massey_decode()
    # test_bch()
    
