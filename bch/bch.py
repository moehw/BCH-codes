import math
import random

from finitefield import get_logarithm_table, get_cyclotomic_cosets, multiply_polynomials, get_polynomial_from_roots, divide_polynomials, check_polynomial_with_root, binary_to_string, is_bit_set, msb, get_positions_of_binary_ones, polynomial_of_argument_to_power

VERBOSE = False
TEST    = False

class BCH(object):
    """Bose–Chaudhuri–Hocquenghem error-correcting code."""

    def __init__(self, n, power, dist, b, primitive_polynomial):
        """
        Constructs a BCH code with the
        specified parameters.

        :param n: length of a code.
        :param power: the power in size of the
        field GF(2ᵖᵒʷᵉʳ).
        :param dist: a constructive distance.
        :param b: a power of element from witch we choose
        cyclotomic cosets for generator polynomial. 
        :param primitive_polynomial: a primitive polynomial
        used to create logarithm table.
        """
        self.dist = dist
        t = (dist - 1) >> 1 # dist = 2*t + 1 - for binary
        self.t = t
        self.n = n
        self.b = b
        self.power = power
        
        print("n: {}\nd: {}\nt: {}\nb: {}\npower: {}".format(n, dist, t, b, power))

        full_cyclotomic_cosets = get_cyclotomic_cosets(power=power)
        if VERBOSE:
            print("All cosets: {}".format("\n".join(binary_to_string(i, False) for i in full_cyclotomic_cosets)))

        beta_powers = [i % (2 ** power - 1) for i in range(b, b + dist - 1)]
        if VERBOSE:
            print("Beta powers: {}".format(beta_powers))

        cyclotomic_cosets = []
        for coset in full_cyclotomic_cosets:
            # find all factors beta**i in this coset if element i is inbetween b and b+2t-1
            for beta_power in beta_powers:
                if is_bit_set(coset, beta_power):
                    cyclotomic_cosets.append(coset)
                    break
        if VERBOSE:
            print("Cosets: {}".format("\n".join(binary_to_string(i, False) for i in cyclotomic_cosets)))

        logarithm_table = get_logarithm_table(power=power, primitive_polynomial=primitive_polynomial)
        generator_polynomial = calculate_generator_polynomial(
            cyclotomic_cosets=cyclotomic_cosets,
            logarithm_table=logarithm_table,
            power=power,
            polynomial_count=len(cyclotomic_cosets)
        )

        if TEST:
            for beta_power in beta_powers:
                generator_with_root = check_polynomial_with_root(generator_polynomial, beta_power, power, logarithm_table)

                if VERBOSE:
                    print("Check generator with root {}: {}".format(beta_power, bin(generator_with_root)[2:]))
                
                assert generator_with_root == 0, \
                    "Root is not root of generator polynomial"

        self.primitive_polynomial = primitive_polynomial
        self.cyclotomic_cosets = cyclotomic_cosets
        self.logarithm_table = logarithm_table
        self.generator_polynomial = generator_polynomial

        k = n - msb(generator_polynomial)
        self.k = k
        print("k: {}".format(k))
        print("Generator polynomial: {} ({:b})".format(binary_to_string(generator_polynomial), generator_polynomial))

    def encode(self, message):
        return encode(bch.generator_polynomial, message)

    def decode(self, message):
        return decode(self.primitive_polynomial, message, self.cyclotomic_cosets, self.logarithm_table, self.power, self.t, self.n, self.k, self.b)


def calculate_generator_polynomial(cyclotomic_cosets, logarithm_table, power, polynomial_count):
    """
    Calculates a generator polynomial of
    a particular BCH code.
    
    :param cyclotomic_cosets: cyclotomic
    cosets of primitive roots of a Galois
    field.
    :param logarithm_table: a convenient form
    of the field elements for multiplication.
    :param power: the power in GF(2ᵖᵒʷᵉʳ).
    :param t: a number of errors to correct.
    
    :returns: a generator polynomial which is
    a product of polynomails obtained from
    cyclotomic cosets.
    """
    generator_polynomial = get_polynomial_from_roots(
                        roots=cyclotomic_cosets[0],
                        power=power,
                        logarithm_table=logarithm_table)
                
    for i in range(1, polynomial_count):
        generator_polynomial = multiply_polynomials(
            polynomial1=generator_polynomial,
            polynomial2=get_polynomial_from_roots(
                roots=cyclotomic_cosets[i],
                power=power,
                logarithm_table=logarithm_table)
        )

    return generator_polynomial

def encode(generator_polynomial, message):
    """
    Encodes a message by shifting
    it to the highest power and
    adding a divided by a generator
    polynomial message to it.
    
    :param generator_polynomial:
    a generator polynomial of a
    certain finite field.
    :param message: a message to encode.
    
    :returns: an encoded message.
    """
    power = msb(generator_polynomial)
    message = message << power
    return message ^ divide_polynomials(polynomial1=message, polynomial2=generator_polynomial)[1]

def decode(primitive_polynomial, received_message, cyclotomic_cosets, logarithm_table, power, t, n, k, b):
    """
    Decodes a received message.
    
    :param primitive_polynomial: a primitive
    polynomial over a Galois field.
    :param received_message: a message is
    to be decoded.
    :param cyclotomic_cosets: cyclotomic
    cosets of primitive roots of a Galois
    field.
    :param logarithm_table: a mapping from
    the power of a primitive element of a Galois
    field to this element.
    :param power: the power in the Galois field
    G(2ᵖᵒʷᵉʳ).
    :param t: a number of errors to be corrected.
    
    :returns: a decoded message.
    """
    syndromes = get_syndromes(
        primitive_polynomial=primitive_polynomial,
        received_message=received_message,
        cyclotomic_cosets=cyclotomic_cosets,
        logarithm_table=logarithm_table,
        power=power,
        t=t,
        b=b)

    sigma = berlekamp_massey_decode(
        syndromes=syndromes,
        logarithm_table=logarithm_table,
        power=power,
        t=t)

    roots = find_roots_of_sigma(
        sigma=sigma,
        power=power,
        logarithm_table=logarithm_table)

    error_positions = get_error_positions(roots=roots, power=power)
    for position in error_positions:
        received_message ^= 1 << position

    received_message >>= (n - k)
    return received_message

def get_syndromes(primitive_polynomial, received_message, cyclotomic_cosets, logarithm_table, power, t, b):
    
    length = t * 2
    syndromes = [0] * length

    beta_powers = [i % (2 ** power - 1) for i in range(b, b + (2 * t))]

    syndromes_t = []
    for i in beta_powers:
        r = check_polynomial_with_root(received_message, 1 << i, 6, logarithm_table)
        syndromes_t.append(r)

    return syndromes_t

def berlekamp_massey_decode(syndromes, logarithm_table, power, t):
    flipped_logarithm_table = flip_dictionary(dictionary=logarithm_table)
    number_of_elements = 2 ** power - 1
    
    capital_lambdas = [[0] * (2 * t + 1)] * (2 * t + 1)
    capital_lambdas[-1][0] = 1
    capital_lambdas[0][0] = 1
    capital_ls = [0] * number_of_elements

    deltas = [0] * (2 * t + 1)
    deltas[0] = 1
    for k in range(1, 2 * t + 1):
        # calculates ∆ (delta)
        for i in range(capital_ls[k - 1] + 1):
            if (k - 1 - i < len(syndromes)) and (syndromes[k - 1 - i] >= 0) and (capital_lambdas[k - 1][i] > 0):
                deltas[k] ^= logarithm_table[
                    (flipped_logarithm_table[capital_lambdas[k - 1][i]] + syndromes[k - 1 - i]) % number_of_elements
                    ]
        # ends
        m_bag = []
        for i in range(k):
            if capital_ls[i] == capital_ls[k - 1]:
                m_bag.append(i)
        m = sorted(m_bag)[0]
        if deltas[k] == 0:
            capital_lambdas[k] = capital_lambdas[k - 1].copy()
            capital_ls[k] = capital_ls[k - 1]
        else:
            capital_lambdas[k] = capital_lambdas[k - 1].copy()
            multiplied_c_l = [0] * (len(capital_lambdas[m - 1]) + k - m)
            for i in range(len(capital_lambdas[m - 1])):
                if capital_lambdas[m - 1][i] > 0:
                    multiplied_c_l[i + k - m] = logarithm_table[
                        (flipped_logarithm_table[capital_lambdas[m - 1][i]] +
                         flipped_logarithm_table[deltas[k]] +
                         number_of_elements - flipped_logarithm_table[deltas[m]]
                         ) % number_of_elements
                        ]
            #                    [k-1]
            # inlines arrays of Λ      and the addition from the k-th iteration,
            #             [k]
            # calculates Λ   .
            who = len(capital_lambdas[k]) - len(multiplied_c_l)
            if who <= 0:
                for i in range(who):
                    capital_lambdas[k][i].append(0)
            else:
                for i in range(who):
                    multiplied_c_l.append(0)
            for i in range(len(capital_lambdas[k])):
                capital_lambdas[k][i] = capital_lambdas[k][i] ^ multiplied_c_l[i]
            # ends
            capital_ls[k] = max(capital_ls[k - 1], capital_ls[m - 1] + k - m)
    
    result = []
    for i in capital_lambdas[2 * t]:
        result.append(flipped_logarithm_table[i])
    return result

def find_roots_of_sigma(sigma, power, logarithm_table):
    """
    Tries all the elements of a field to solve an
    error polynomial equls to zero.
    
    :param sigma: a sigma is an array representing
    an error locator polynomial.
                               2        14    2
    The result looks like 1 + α  * x + α   * x .
    The array represents the powers of x accordingly
                                  0               1
    to the position. array[0] is x , array[1] is x .
    In these cells the powers of α are stored.
    :param power: the power in a Galois field GF(2ᵖᵒʷᵉʳ).
    :param logarithm_table: a mapping of a power
    to a corresponding element of a particular
    Galois field.
    
    :returns: an array of roots of a sigma.
    """
    roots = []
    for candidate in range(2 ** power - 1):
        result = logarithm_table[sigma[0]]
        for polynomial_power in range(1, get_order_of_sigma(sigma=sigma) + 1):
            if sigma[polynomial_power] >= 0:
                result ^= logarithm_table[(sigma[polynomial_power] + candidate * polynomial_power) % (2 ** power - 1)]
        if result == 0:
            roots.append(candidate)
    return roots

def get_error_positions(roots, power):
    """
    Flips the established roots of
    an error locator polynomial to
    the positions of errors.
    
    :param roots: an array of roots
    of an error locator polynomial.
    :param power: the power in a Galois
    field GF(2ᵖᵒʷᵉʳ).
    
    :returns: an array of error
    positions.
    """
    positions = []
    for root in roots:
        positions.append((2 ** power - 1 - root) % (2 ** power - 1))
    return positions

def get_order_of_sigma(sigma):
    """
    Returns an order of an
    arror locator polynomial.
    
    :param sigma: an error
    locator polynomial.
    
    :returns: an order of
    sigma as an integer.
    """
    counter = 0
    for i in reversed(sigma):
        if i > -1:
            return len(sigma) - 1 - counter
        counter += 1
    return 0

def flip_dictionary(dictionary):
    """
    Flips a dictionary: key <-> value.
    
    :param dictionary: a dictionary
    to flip values and keys.
    
    :returns: a flipped dictionary.
    """
    return dict((v, k) for k, v in dictionary.items())

# Message input

def message_to_bits_and_split_on_blocks(message, n):
    """
    Translates a message to bit representation and splits bits
    on blocks of length n.
    
    :param message: a message to translate.
    :param n: length of a block.
    
    :returns: an array of binary blocks.
    """
    bits = text_to_bits(message)
    blocks = []

    while bits != 0:
        block = bits & (2 ** n - 1)
        blocks.append(block)
        bits >>= n
    return blocks

def text_to_bits(text, encoding='utf-8', errors='surrogatepass'):
    """
    Translates given text to bits.
    
    :param text: a string to translate.
    :param encoding: a used encoding.
    :param errors: an error attitude.
    
    :returns: a binary vector.
    """
    return int.from_bytes(text.encode(encoding, errors), 'big')

def bits_to_message_and_concat_blocks(blocks, n):
    """
    Translates an array of bit vectors to string.
    
    :param blocks: an array of bit vectors.
    :param n: length of a bit vector.
    
    :returns: a string of united translated vectors.
    """
    result = 0
    for block in reversed(blocks):
        result <<= n
        result ^= block
    return result, text_from_bits(result)

def text_from_bits(bits, encoding='utf-8', errors='surrogatepass'):
    """
    Translates a bit array to string.
    
    :param bits: a bit array.
    :param encoding: a used encoding.
    :param errors: an error attitude.
    
    :returns: a string.
    """
    return bits.to_bytes((bits.bit_length() + 7) // 8, 'big').decode(encoding, errors) or '\0'

# Samples

def lab():
    primpoly = 0b1101101
    n = 63
    d = 11
    t = (d - 1) >> 1
    b = 55
    power = 6
    bch = BCH(n, power, d, b, primpoly)
    return bch

def lecture():
    bch = BCH(15, 4, 5, 1, 0b10011)
    return bch

def ksu():
    primpoly = 0b1110011
    n = 63
    d = 8
    t = (d - 1) >> 1
    b = 0
    a = 0b111011111

    bch = BCH(n, 6, d, b, primpoly)
    return bch

def ilya():
    primpoly = 0b1110011
    n = 63
    d = 11
    t = (d - 1) >> 1
    b = 1
    a = 0b10100010010

    print("Input: {} ({:b})".format(binary_to_string(a), a))

    bch = BCH(n, 6, d, b, primpoly)
    res = encode(bch.generator_polynomial, a)
    print("Result: {} ({:b})".format(binary_to_string(res), res))

    res = decode(primpoly, res, bch.cyclotomic_cosets, bch.logarithm_table, 6, t)
    print("Result: {} ({:b})".format(binary_to_string(res), res))

if __name__ == "__main__":
    
    message = 0b11010000101
    bch = lab()

    # message = 0b1010101
    # bch = lecture()

    print("Input: {} ({:b})".format(binary_to_string(message), message))

    result = bch.encode(message)
    print("Result encode: {} ({:b})".format(binary_to_string(result), result))

    # result ^= 0b10000000000000000000000000000000000000
    # print(bin(result))

    result = bch.decode(result)
    print("Result decode: {} ({:b})".format(binary_to_string(result), result))

    assert message == result, \
        "Message was encoded or decoded incorrect"
    