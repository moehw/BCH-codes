from itertools import combinations

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
    assert msb(primitive_polynomial) == power, \
     "The primitive polynomial {:b} is not of the specified power n = {}".format(primitive_polynomial, power)

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
