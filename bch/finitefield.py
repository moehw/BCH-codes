from itertools import combinations

def get_primitive_polynomial(power):
    """
    Retrieves a table of primitive
    polynomials from the given in
    the config file, and returns a
    corresponding to the
    parameters polynomial.
    
    :param power: the power of 2 in
    the length of a code (2ᵖᵒʷᵉʳ).
    
    :returns: a primitive polynomial
    in a binary representation.
    """
    
    prim_polynomials = {
    2 : 0b111,
    3 : 0b1011,
    4 : 0b10011,
    5 : 0b100101,
    6 : 0b1000011,
    7 : 0b10001001,
    8 : 0b100011101,
    9 : 0b1000010001,
    10 : 0b10000001001,
    11 : 0b100000000101,
    12 : 0b1000001010011,
    13 : 0b10000000011011,
    14 : 0b100010001000011,
    15 : 0b1000000000000011,
    16 : 0b10001000000001011,
    17 : 0b100000000000001001,
    18 : 0b1000000000010000001,
    19 : 0b10000000000000100111,
    20 : 0b100000000000000001001,
    # 21 : 0b x^21+x^2+1,
    # 22 : 0b x^22+x+1,
    # 23 : 0b x^23+x^5+1,
    # 24 : 0b x^24+x^7+x^2+x+1,
    # 25 : 0b x^25+x^3+1,
    # 26 : 0b x^26+x^6+x^2+x+1,
    # 27 : 0b x^27+x^5+x^2+x+1,
    # 28 : 0b x^28+x^3+1,
    # 29 : 0b x^29+x^2+1,
    # 30 : 0b x^30+x^23+x^2+x+1,
    # 31 : 0b x^31+x^3+1,
    # 32 : x^32+x^22+x^2+x+1,
    # 33 : x^33+x^13+1,
    # 34 : x^34+x^8+x^4+x^3+1,
    # 35 : x^35+x^2+1,
    # 36 : x^36+x^11+1,
    # 37 : x^37+x^6+x^4+x+1,
    # 38 : x^38+x^6+x^5+x+1,
    # 39 : x^39+x^4+1,
    # 40 : x^40+x^5+x^4+x^3+1,
    # 41 : x^41+x^3+1,
    # 42 : x^42+x^7+x^4+x^3+1,
    # 43 : x^43+x^6+x^4+x^3+1,
    # 44 : x^44+x^6+x^5+x^2+1,
    # 45 : x^45+x^4+x^3+x+1,
    # 46 : x^46+x^8+x^7+x^6+1,
    # 47 : x^47+x^5+1,
    # 48 : x^48+x^9+x^7+x^4+1,
    # 49 : x^49+x^9+1,
    # 50 : x^50+x^4+x^3+x^2+1,
    # 51 : x^51+x^6+x^3+x+1,
    # 52 : x^52+x^3+1,
    # 53 : x^53+x^6+x^2+x+1,
    # 54 : x^54+x^8+x^6+x^3+1,
    # 55 : x^55+x^24+1,
    # 56 : x^56+x^7+x^4+x^2+1,
    # 57 : x^57+x^7+1,
    # 58 : x^58+x^19+1,
    # 59 : x^59+x^7+x^4+x^2+1,
    # 60 : x^60+x+1,
    # 61 : x^61+x^5+x^2+x+1,
    # 62 : x^62+x^6+x^5+x^3+1,
    # 63 : x^63+x+1,
    # 64 : x^64+x^4+x^3+x+1,
    # 65 : x^65+x^18+1,
    # 66 : x^66+x^9+x^8+x^6+1,
    # 67 : x^67+x^5+x^2+x+1,
    # 68 : x^68+x^9+1,
    # 69 : x^69+x^6+x^5+x^2+1,
    # 70 : x^70+x^5+x^3+x+1,
    # 71 : x^71+x^6+1,
    # 72 : x^72+x^10+x^9+x^3+1,
    # 73 : x^73+x^25+1,
    # 74 : x^74+x^7+x^4+x^3+1,
    # 75 : x^75+x^6+x^3+x+1,
    # 76 : x^76+x^5+x^4+x^2+1,
    # 77 : x^77+x^6+x^5+x^2+1,
    # 78 : x^78+x^7+x^2+x+1,
    # 79 : x^79+x^9+1,
    # 80 : x^80+x^9+x^4+x^2+1,
    # 81 : x^81+x^4+1,
    # 82 : x^82+x^9+x^6+x^4+1,
    # 83 : x^83+x^7+x^4+x^2+1,
    # 84 : x^84+x^13+1,
    # 85 : x^85+x^8+x^2+x+1,
    # 86 : x^86+x^6+x^5+x^2+1,
    # 87 : x^87+x^13+1,
    # 88 : x^88+x^11+x^9+x^8+,
    # 89 : x^89+x^38+1,
    # 90 : x^90+x^5+x^3+x^2+1,
    # 91 : x^91+x^8+x^5+x+1,
    # 92 : x^92+x^6+x^5+x^2+1,
    # 93 : x^93+x^2+1,
    # 94 : x^94+x^21+1,
    # 95 : x^95+x^11+1,
    # 96 : x^96+x^10+x^9+x^6+1,
    # 97 : x^97+x^6+1,
    # 98 : x^98+x^11+1,
    # 99 : x^99+x^7+x^5+x^4+1,
    # 100 : x^100+x^37+1,
    # 101 : x^101+x^7+x^6+x+1,
    # 102 : x^102+x^6+x^5+x^3+1,
    # 103 : x^103+x^9+1,
    # 104 : x^104+x^11+x^10+x+1,
    # 105 : x^105+x^16+1,
    # 106 : x^106+x^15+1,
    # 107 : x^107+x^9+x^7+x^4+1,
    # 108 : x^108+x^31+1,
    # 109 : x^109+x^5+x^4+x^2+1,
    # 110 : x^110+x^6+x^4+x+1,
    # 111 : x^111+x^10+1
    }

    return prim_polynomials[power]

def get_logarithm_table(power, primitive_polynomial):
    """
    Builds a logarithmic table where a
    logarithm is mapped to the binary
    representation of the corresponding
    polynomial.
    
    :param primitive_polynomial.
    :param power: the power in size of the
    field GF(2ᵖᵒʷᵉʳ).

    :returns: a logarithm table of the field.
    """
    if msb(primitive_polynomial) != power:
        raise ValueError("The primitive polynomial {:b} is not of the specified power n = {}".format(primitive_polynomial, power))

    logarithm_table = {-1: 0}

    for i in range(power):
        logarithm_table[i] = 1 << i

    logarithm_table[power] = trim_polynomial(polynomial=primitive_polynomial, length=power)
    
    for i in range(power + 1, 2 ** power - 1):
        multiplied_by_x_polynomial = logarithm_table[i - 1] << 1
        if multiplied_by_x_polynomial & (2 ** power):
            multiplied_by_x_polynomial ^= logarithm_table[power]
        logarithm_table[i] = trim_polynomial(polynomial=multiplied_by_x_polynomial, length=power)
    return logarithm_table

def multiply_polynomials(polynomial1, polynomial2):
    """
    Multiplies two polynomials in GF(2ᵖᵒʷᵉʳ).

    :param polynomial1: 1st polynomial.
    :param polynomial2: 2nd polynomial.
    
    :returns: the product of two polynomials.
    """
    result = 0
    for i in range(msb(polynomial2) + 1):
        if is_bit_set(polynomial2, i):
            result ^= polynomial1 << i
    return result

def divide_polynomials(polynomial1, polynomial2):
    """
        101010100000000 | 111010001
       _________________|________________
                . . .   |   1110101 <- quotient
              __________|
               11100101 <- reminder

    :param polynomial1: 1st polynomial.
    :param polynomial2: 2nd polynomial.

    :returns: the quotient and the remainder.
    """
    quotient = 0
    reminder = polynomial1
    while msb(reminder) >= msb(polynomial2):
        shift = msb(reminder) - msb(polynomial2)
        reminder ^= polynomial2 << shift
        quotient ^= 1 << shift
    return quotient, reminder

def trim_polynomial(polynomial, length):
    """
    Cuts off extra bits form the
    polynomial.

    :param polynomial: a polynomial
    to trim.
    :param length: length of the
    polynomial.
    
    :returns: a reduced polynomial.
    """
    return polynomial & ((2 ** length) - 1)

def polynomial_of_argument_to_power(polynomial, power):
    """
    Calculates terms of a polynomial putting
    in corresponding power of the variable
    and reduces the resulting polynomial.
    
    :param polynomial: a polynomial to
    calculate.
    :param power: the power to put in.
    
    :returns: a reduced polynomial of the
    specified power.
    """
    
    length = msb(polynomial) + 1
    result = 0

    for i in range(length):
        if is_bit_set(polynomial, i):
            result |= 1 << i * power
    
    return result

def get_cyclotomic_cosets(power):
    """
    Fills a list of cyclotomic cosets.
    For GF(16) = GF(2⁴) cyclotomic
    cosets are:
    m0(x)
    m1(x) = m2(x) = m4(x) = m8(x) 
    m3(x) = m6(x) = m9(x) = m12(x)
    m5(x) = m10(x)
    m7(x) = m11(x) = m13(x) = m14(x)

    :param power: the power in size of the
    field GF(2ᵖᵒʷᵉʳ).

    :returns: a list of cyclotomic cosets.
    """
    cyclotomic_cosets = []
    all_cyclotomic_members = 1
    i = 0

    while all_cyclotomic_members < 2 ** (2 ** power - 1) - 1:
        cyclotomic_cosets.append(0)
        k = 0
        while True:
            if not is_bit_set(all_cyclotomic_members, k):
                break
            k += 1
        while True:
            k = k % (2 ** power - 1)
            if is_bit_set(cyclotomic_cosets[i], k):
                break
            cyclotomic_cosets[i] ^= 1 << k
            k *= 2
        all_cyclotomic_members ^= cyclotomic_cosets[i]
        i += 1

    cyclotomic_cosets.append(1) # {a0}
    
    return cyclotomic_cosets

def check_polynomial_with_root(polynomial, root, logarithm_table):
    """
    Substitutes roots in polynomial.

    :param polynomial: a polynomial to
    calculate.
    :param root: a binary root,
    where positions of 1s mean
    the power of the field.
    :param logarithm_table:
    a table which maps logarithms
    to polynomials - members of
    the field.
    
    :returns: a binary vector.
    """
    number_of_field_elements = len(logarithm_table) - 1
    power_array = get_positions_of_binary_ones(polynomial)

    result = 0
    for i in power_array:
        result ^= logarithm_table[(i * root) % number_of_field_elements]
    
    return result

def get_polynomial_from_roots(roots, logarithm_table):
    """
    Performs multiplication of a
    polynomial represented by its
    roots in form:
    (x - aᵏ¹)*(x - aᵏ²)...*(x - aᵏⁿ).
    
    :param roots: a binary vector
    of roots, where positions of
    1s mean the power a primitive
    element a of the field.
    :param logarithm_table:
    a table which maps logarithms
    to polynomials - members of
    the field.
    
    :returns: a binary vector
    represents a polynomial in
    form of xˡ¹ + xˡ² ... + xˡⁿ + 1
    """
    if roots == 0:
        return 0

    number_of_field_elements = len(logarithm_table) - 1
    root_array = get_positions_of_binary_ones(number=roots)
    polynomial = 1 << len(root_array)
    for i in range(len(root_array)):
        coefficient = 0
        for combination in combinations(root_array, i + 1):
            coefficient ^= logarithm_table[sum(combination) % number_of_field_elements]
        addition = coefficient << len(root_array) - i - 1
        polynomial ^= addition

    return polynomial

def get_positions_of_binary_ones(number):
    """
    Gets positions of 1s in a binary
    vector.

    :param number: a vector to
    analyze.

    :returns: positions of 1s.
    """
    length = msb(number) + 1
    result = []

    for i in range(0, length):
        if is_bit_set(number, i):
            result.append(i)

    return result

def is_bit_set(polynomial, n):
    """
    Checks if bit is setted.

    :param polynomial: a polynomial
    to check.
    :param n: a position of bit to
    check.

    :returns: a status of bit.
    """
    return (polynomial & (1 << n))

def binary_to_string(polynomial, direction_polynomial=True):
    """
    Represents binary polynomial in
    user-friendly form.

    :param polynomial: a polynomial to
    convert.
    :param direction_polynomial: if true,
    use "+" delimeter and starts from higher
    power. Otherwise, delimeter is "," and
    from lower power. 

    :returns: a string of polynomial.
    """
    
    if polynomial == 0:
        return "0"

    superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    s = "" if not is_bit_set(polynomial, 0) else "1"
    for i in range(1, msb(polynomial) + 1):
        if is_bit_set(polynomial, i):
            if direction_polynomial:
                s = "x" + ("" if i == 1 else str(i).translate(superscript)) + \
                    ("" if s == "" else " + ") + s
            else:
                s += ("" if s == "" else ", ") + \
                    "x" + ("" if i == 1 else str(i).translate(superscript))
    return s

def msb(number):
    """
    Most significant bit (MSB).
    Position starts from zero.
    
    :param number: a binary
    number to check.

    :returns: a position of 
    most significant bit.
    """

    position = 0
    number = number >> 1
    
    while number != 0:
        position += 1
        number = number >> 1
    
    return position
