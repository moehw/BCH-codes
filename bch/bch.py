import math
import random

from finitefield import *

VERBOSE = True
TEST    = False

class BCH:
    """Bose–Chaudhuri–Hocquenghem error-correcting code."""

    def __init__(self, n, dist, b, primitive_polynomial):
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
        power = msb(primitive_polynomial)
        self.power = power
        
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
                generator_with_root = check_polynomial_with_root(generator_polynomial, beta_power, logarithm_table)

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

        print("n: {}, d: {}, t: {}, b: {}, power: {}, k: {}".format(n, dist, t, b, power, k))
        print("Generator polynomial: {} ({:b})".format(binary_to_string(generator_polynomial), generator_polynomial))

    def encode(self, message):
        return encode(self.generator_polynomial, message)

    def decode(self, message):
        return decode(self.primitive_polynomial, message, self.cyclotomic_cosets, self.logarithm_table, self.power, self.t, self.n, self.k, self.b)

def create_code_with_fix_speed(block_size, speed, dist):
    """
    :param block_size: a size of one block of
    input data to be encoded/decoded.
    :param speed: a speed of code, equal to
    k divided by n.
    :param dist: a disigned distination between
    codewords.

    :returns: a BCH code class with given
    code speed.
    """

    k = round(block_size * speed)

    return BCH(block_size, dist, 1, get_primitive_polynomial(block_size - k - dist + 1, 1))

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
                        logarithm_table=logarithm_table)
                
    for i in range(1, polynomial_count):
        generator_polynomial = multiply_polynomials(
            polynomial1=generator_polynomial,
            polynomial2=get_polynomial_from_roots(
                roots=cyclotomic_cosets[i],
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
    message <<= power
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
    
    :returns: a decoded message and status.
    """
    syndromes, is_error = get_syndromes(
        primitive_polynomial=primitive_polynomial,
        received_message=received_message,
        cyclotomic_cosets=cyclotomic_cosets,
        logarithm_table=logarithm_table,
        power=power,
        t=t,
        b=b)
    
    if not is_error:
        if VERBOSE:
            print("Message has no errors")
        return received_message

    sigma = berlekamp_massey_decode(
        syndromes=syndromes,
        logarithm_table=logarithm_table,
        power=power,
        t=t,
        primitive_polynomial=primitive_polynomial)

    roots = find_roots_of_sigma(
        sigma=sigma,
        power=power,
        logarithm_table=logarithm_table)

    error_positions = get_error_positions(roots=roots, power=power)
    for position in error_positions:
        if VERBOSE:
            print("Error in position {}".format(position))
        received_message ^= 1 << position

    received_message >>= (n - k)
    return received_message

def get_syndromes(primitive_polynomial, received_message, cyclotomic_cosets, logarithm_table, power, t, b):
    """
    Calculates syndromes as values of received
    message at zeroes of generator polynomial.

    Calculates syndromes based on a particular
    received message, in number of 2 * t.
    E.g. received message is 100100100000001
    should be divided by each of the minimal
    polynomials. The reminder is 1110 if it is
    divided by M1 = 11001. Using roots of M1
    a, a² , a⁴ three syndromes are obtained:
    S1 = 1110, due to x = a;
                                
    S2 = 101010100, due to x = a² for 1110;
    S3 = 10001000100010000, due to x = a³ for
    1110.
    Note that 1110 means x³ + x² + x¹.

    :param primitive_polynomial: a primitive
    polynomial, primitive in sense of the
    given GF(2ᵖᵒʷᵉʳ).
    :param received_message: a message which
    was received by the decoder.
    :param cyclotomic_cosets: cyclotomic cosets
    for building minimal polynomials.
    :param logarithm_table: a convenient form
    of the field elements for multiplication.
    :param power: the power in GF(2ᵖᵒʷᵉʳ).
    :param t: a number of errors to correct.
    
    :returns: a list of powers of a primitive
    element a as shortcuts for polynomials.
    """
    length = t * 2
    # syndromes = [0] * length
    syndromes = []
    is_error = False

    beta_powers = [i % (2 ** power - 1) for i in range(b, b + (2 * t))]

    flipped_logarithm_table = flip_dictionary(dictionary=logarithm_table)
    
    for i in beta_powers:
        syndrome_polynomial = received_message

        for coset in cyclotomic_cosets:
            if is_bit_set(coset, i):
                syndrome_polynomial = divide_polynomials(
                    polynomial1=syndrome_polynomial,
                    polynomial2=get_polynomial_from_roots(
                        roots=coset,
                        logarithm_table=logarithm_table)
                )[1]
                break

        syndrome_polynomial_of_argument_to_power = polynomial_of_argument_to_power(
            polynomial=syndrome_polynomial,
            power=i)

        if msb(syndrome_polynomial_of_argument_to_power) >= msb(primitive_polynomial):
            syndrome_polynomial_of_argument_to_power = divide_polynomials(
                polynomial1=syndrome_polynomial_of_argument_to_power,
                polynomial2=primitive_polynomial
            )[1]

        syndrome = flipped_logarithm_table[syndrome_polynomial_of_argument_to_power]

        if syndrome != -1:
            is_error = True

        syndromes.append(syndrome)
        
    return syndromes, is_error

def berlekamp_massey_decode(syndromes, logarithm_table, power, t, primitive_polynomial):
    """
    Calculates an error locator polynomial using
    the Berlekamp-Massey algorithm.
    :param syndromes: a calculated array of
    syndromes corresponding to the received
    message with error vector.
    :param logarithm_table: a mapping of a
    power to a corresponding element of a particular
    Galois field.
    :param power: the power in a Galois field GF(2ᵖᵒʷᵉʳ).
    :param t: a number of error to be corrected.
    
    :returns: an array representing an error locator
    polynomial.
                               2        14    2
    The result looks like 1 + α  * x + α   * x .
    The array represents the powers of x accordingly
                                  0               1
    to the position. array[0] is x , array[1] is x .
    In these cells the powers of α are stored.
    """
    
    flipped_logarithm_table = flip_dictionary(dictionary=logarithm_table)
    number_of_elements = 2 ** power - 1
    
    C = [0] * (2 * t) # error locator polynomial
    C[0] = 1
    L = 0 # len of C

    B = [0] * (2 * t) # copy of C on last L update
    B[0] = 1

    for n in range(2 * t):

        # ∆ (delta)
        delta = 0
        for i in range(0, L + 1):
            if not (C[i] and logarithm_table[syndromes[n - i]]):
                # one of the factors is zero
                continue

            C_idx = flipped_logarithm_table[C[i]]
            s_idx = syndromes[n - i]
            mult = (C_idx + s_idx) % number_of_elements
            delta ^= logarithm_table[mult] 

        if VERBOSE:
            print("\nn: {}, delta: {}".format(n, delta))

        for i in range(len(B) - 1, -1, -1):
            B[i] = 0 if (i == 0) else B[i - 1] # x*B(x)
        
        if delta == 0:
            # L not updates
            pass
        else:
            T = C.copy() # temporary polinomial

            for i in range(len(C)):
                b_idx = flipped_logarithm_table[B[i]]
                delta_idx = flipped_logarithm_table[delta]
                if b_idx != -1:
                    mult = (b_idx + delta_idx) % number_of_elements
                    T[i] ^= logarithm_table[mult]

            if 2 * L <= n:
                # L updates, nullify delta

                for i in range(len(C)):
                    B[i] = 0
                    delta_idx = number_of_elements - flipped_logarithm_table[delta]
                    C_idx = flipped_logarithm_table[C[i]]
                    if C_idx != -1:
                        div = (C_idx + delta_idx) % number_of_elements
                        B[i] = logarithm_table[div]

                C = T.copy()

                L = n + 1 - L
            else:
                C = T.copy()

    result = []
    for poly in C:
        result.append(flipped_logarithm_table[poly])
    return result

def find_roots_of_sigma(sigma, power, logarithm_table):
    """
    Tries all the elements of a field to solve an
    error polynomial equls to zero. Chien algorithm.
    
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
    error locator polynomial.
    
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

if __name__ == '__main__':
    bch = BCH(31, 2 * 3 + 1, 1, get_primitive_polynomial(6, 1))
    block = 5211
    encode_ = 1366166387
    distorted = 1131154291

    exp_encode = bch.encode(block)

    assert exp_encode == encode_

    exp_block = bch.decode(distorted)

    assert exp_block == block