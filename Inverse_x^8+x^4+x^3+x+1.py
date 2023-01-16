import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control

def Inverison_Itoh_Tsujii_based(eng) :

    n = 8
    a = eng.allocate_qureg(n)
    a1 = eng.allocate_qureg(n)

    check = eng.allocate_qureg(n)
    if (resource_check != 1):
        Round_constant_XOR(eng, a, 0xab, n)
        Round_constant_XOR(eng, check, 0xab, n)

    ancilla  = eng.allocate_qureg(38)
    count = 0

    copy(eng, a, a1, n) # a1 = a
    a1 = Squaring(eng, a1, n) # a1 = a^2

    a2 = []
    a2, count, ancilla = recursive_karatsuba(eng, a, a1, n, count, ancilla)
    a2 = Reduction(eng, a2)

    # (a * a^2) ^64
    a = Squaring(eng, a, n)
    a = Squaring(eng, a, n)
    a = Squaring(eng, a, n)
    a = Squaring(eng, a, n)
    a = Squaring(eng, a, n)
    a = Squaring(eng, a, n) #^64

    count = 0
    a3 = []
    a3, count, ancilla = recursive_karatsuba(eng, a2, a, n, count, ancilla)
    a3 = Reduction(eng, a3)

    a2 = Squaring(eng, a2, n)
    a2 = Squaring(eng, a2, n)

    count = 0
    a4 = []
    a4, count, ancilla = recursive_karatsuba(eng, a2, a3, n, count, ancilla)
    a4 = Reduction(eng, a4)

    a2 = Squaring(eng, a2, n)
    a2 = Squaring(eng, a2, n)

    count = 0
    a5 = []
    if (resource_check != 1):
        a5, count, ancilla = recursive_karatsuba(eng, a2, a4, n, count, ancilla) #00011100
    else:
        a5, count, ancilla = recursive_karatsuba_omit_reverse(eng, a2, a4, n, count, ancilla)  # 00011100
    a5 = Reduction(eng, a5)

    a5 = Squaring(eng, a5, n)

    if(resource_check != 1):
        print_state(eng, a5, n)

    if (resource_check != 1):
        count = 0
        result = []
        result, count, ancilla = recursive_karatsuba(eng, a5, check, n, count, ancilla)
        result = Reduction(eng, result)
        print('a * a^-1')
        print_state(eng, result, n)

def copy(eng, a, b, n):
    for i in range(n):
        CNOT | (a[i], b[i])

def Reduction(eng, result):
    CNOT | (result[8], result[0])
    CNOT | (result[12], result[0])
    CNOT | (result[13], result[0])

    CNOT | (result[12], result[1])
    CNOT | (result[13], result[2])
    CNOT | (result[8], result[3])

    CNOT | (result[11], result[4])
    CNOT | (result[12], result[5])
    CNOT | (result[11], result[7])

    #######
    CNOT | (result[14], result[8])  # 8 = green
    CNOT | (result[9], result[8])
    CNOT | (result[10], result[9])  # 9 = yellow

    CNOT | (result[12], result[14])  # 14 = blue
    CNOT | (result[13], result[10])  # 10 = orange
    CNOT | (result[11], result[10])

    #######
    CNOT | (result[14], result[3])
    CNOT | (result[14], result[7])
    CNOT | (result[10], result[3])

    CNOT | (result[10], result[6])
    CNOT | (result[8], result[1])
    CNOT | (result[8], result[4])

    CNOT | (result[9], result[2])
    CNOT | (result[9], result[5])

    return result[0:8]

def Squaring(eng, vector, n):
    P = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    L = [1, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 0, # (x4, x6) (x5, x6)
         0, 0, 0, 0, 0, 0, 1, 1] # (x6, x7)

    U = [1, 0, 0, 0, 1, 0, 1, 0, # (x4, x0) (x6, x0)
         0, 1, 0, 0, 0, 1, 0, 0, # (x5, x1)
         0, 0, 1, 0, 1, 0, 0, 1, # (x4, x2) (x7, x2)
         0, 0, 0, 1, 0, 1, 0, 0, # (x5, x3)
         0, 0, 0, 0, 1, 0, 1, 1, # (x6, x4) (x7, x4)
         0, 0, 0, 0, 0, 1, 1, 0, # (x6, x5)
         0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 0, 0, 0, 0, 0, 1]

    Apply_U(eng, vector, n, U)
    Apply_L(eng, vector, n, L)

    out = []
    out.append(vector[0])
    out.append(vector[4])
    out.append(vector[1])
    out.append(vector[6])

    out.append(vector[2])
    out.append(vector[5])
    out.append(vector[3])
    out.append(vector[7])

    return out


def Apply_U(eng, key, n, U):
    for i in range(n - 1):
        for j in range(n - 1 - i):
            if (U[(i * n) + 1 + i + j] == 1):
                CNOT | (key[1 + i + j], key[i])


def Apply_L(eng, key, n, L):
    for i in range(n - 1):
        for j in range(n - 1, 0 + i, -1):
            if (L[n * (n - 1 - i) + (n - 1 - j)] == 1):
                CNOT | (key[n - 1 - j], key[n - 1 - i])


def recursive_karatsuba(eng, a, b, n, count, ancilla): #n=4

    if(n==1):
        c = eng.allocate_qubit()
        Toffoli_gate(eng, a, b, c)

        return c, count, ancilla

    c_len = 3**math.log(n, 2) #9 #3
    r_low = n//2    #2 #1

    if(n%2!=0):
        r_low = r_low +1 # n=3 -> 2, n=4 -> 2

    r_a = []
    r_b = []

    # Provide rooms and prepare operands
    r_a = ancilla[count:count + r_low]

    count = count + r_low

    r_b = ancilla[count:count + r_low]

    count = count + r_low

    with Compute(eng):
        for i in range(r_low):
            CNOT | (a[i], r_a[i])
        for i in range(n//2):
            CNOT | (a[r_low + i], r_a[i])
        for i in range(r_low):
            CNOT | (b[i], r_b[i])
        for i in range(n//2):
            CNOT | (b[r_low + i], r_b[i])

    # upper-part setting
    if(r_low == 1):
        c = eng.allocate_qureg(3)
        Toffoli_gate(eng, a[0], b[0], c[0])
        Toffoli_gate(eng, a[1], b[1], c[2])
        CNOT | (c[0], c[1])
        CNOT | (c[2], c[1])
        Toffoli_gate(eng, r_a, r_b, c[1])

        Uncompute(eng)
        return c, count, ancilla

    c_a = []
    c_b = []
    c_r = []

    c_a, count, ancilla = recursive_karatsuba(eng, a[0:r_low], b[0:r_low], r_low, count, ancilla)# 2 qubits     # 0~2
    c_b, count, ancilla = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n//2, count, ancilla)#2 qubits        # 3~5
    c_r, count, ancilla = recursive_karatsuba(eng, r_a[0:r_low], r_b[0:r_low], r_low, count, ancilla) #2qubits  # 6~8

    Uncompute(eng)

    result = []
    result = combine(eng, c_a, c_b, c_r, n)

    return result, count, ancilla

def recursive_karatsuba_omit_reverse(eng, a, b, n, count, ancilla): #n=4

    if(n==1):
        c = eng.allocate_qubit()
        Toffoli_gate(eng, a, b, c)

        return c, count, ancilla

    c_len = 3**math.log(n, 2) #9 #3
    r_low = n//2    #2 #1

    if(n%2!=0):
        r_low = r_low +1 # n=3 -> 2, n=4 -> 2

    r_a = []
    r_b = []

    # Provide rooms and prepare operands
    r_a = ancilla[count:count + r_low]

    count = count + r_low

    r_b = ancilla[count:count + r_low]

    count = count + r_low

    with Compute(eng):
        for i in range(r_low):
            CNOT | (a[i], r_a[i])
        for i in range(n//2):
            CNOT | (a[r_low + i], r_a[i])
        for i in range(r_low):
            CNOT | (b[i], r_b[i])
        for i in range(n//2):
            CNOT | (b[r_low + i], r_b[i])

    # upper-part setting
    if(r_low == 1):
        c = eng.allocate_qureg(3)
        Toffoli_gate(eng, a[0], b[0], c[0])
        Toffoli_gate(eng, a[1], b[1], c[2])
        CNOT | (c[0], c[1])
        CNOT | (c[2], c[1])
        Toffoli_gate(eng, r_a, r_b, c[1])


        return c, count, ancilla

    c_a = []
    c_b = []
    c_r = []

    c_a, count, ancilla = recursive_karatsuba(eng, a[0:r_low], b[0:r_low], r_low, count, ancilla)# 2 qubits     # 0~2
    c_b, count, ancilla = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n//2, count, ancilla)#2 qubits        # 3~5
    c_r, count, ancilla = recursive_karatsuba(eng, r_a[0:r_low], r_b[0:r_low], r_low, count, ancilla) #2qubits  # 6~8

    result = []
    result = combine(eng, c_a, c_b, c_r, n)

    return result, count, ancilla

def combine(eng, a, b, r, n):
    if (n % 2 != 0):
        # n = 13########
        for i in range(n):
            CNOT | (a[i], r[i])
        for i in range(n - 2):
            CNOT | (b[i], r[i])

        for i in range(n // 2):
            CNOT | (a[n // 2 + 1 + i], r[i])
        for i in range(n // 2):
            CNOT | (b[i], r[n // 2 + 1 + i])

        out = []
        for i in range(n // 2 + 1):  # (2n-1) = n//2 + 1 + n ? / 13 = 3+1+7+?
            out.append(a[i])
        for i in range(n):
            out.append(r[i])
        for i in range((2 * n - 1) - n // 2 - 1 - n):
            out.append(b[n // 2 + i])

        return out

    half_n = int(n/2) #n=4
    for i in range(n-1):
        CNOT | (a[i], r[i])
        CNOT | (b[i], r[i])
    for i in range(half_n-1):
        CNOT | (a[half_n+i], r[i])
        CNOT | (b[i], r[half_n+i])

    result = []
    for i in range(half_n):
        result.append(a[i])
    for i in range(n-1):
        result.append(r[i])
    for i in range(half_n):
        result.append(b[half_n-1+i])

    return result

def room(eng, length):

    room = eng.allocate_qureg(length)

    return room

def copy(eng, a, b, length):
    for i in range(length):
        CNOT | (a[i], b[i])

def Toffoli_gate(eng, a, b, c):

    if(NCT):
        Toffoli | (a, b, c)
    else:
        if (resource_check):
            if(AND_check):
                ancilla = eng.allocate_qubit()
                H | c
                CNOT | (b, ancilla)
                CNOT | (c, a)
                CNOT | (c, b)
                CNOT | (a, ancilla)
                Tdag | a
                Tdag | b
                T | c
                T | ancilla
                CNOT | (a, ancilla)
                CNOT | (c, b)
                CNOT | (c, a)
                CNOT | (b, ancilla)
                H | c
                S | c

            else:
                Tdag | a
                Tdag | b
                H | c
                CNOT | (c, a)
                T | a
                CNOT | (b, c)
                CNOT | (b, a)
                T | c
                Tdag | a
                CNOT | (b, c)
                CNOT | (c, a)
                T | a
                Tdag | c
                CNOT | (b, a)
                H | c


def Round_constant_XOR(eng, k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X | k[i]

def print_state(eng, b, n):
    All(Measure) | b
    print('Result : ', end='')
    for i in range(n):
        print(int(b[n-1-i]), end='')
    print('\n')

global resource_check
global AND_check
global NCT

NCT = 1
resource_check = 0
classic = ClassicalSimulator()
eng = MainEngine(classic)
Inverison_Itoh_Tsujii_based(eng)
eng.flush()
print('\n')

resource_check = 1
NCT = 0
AND_check = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
Inverison_Itoh_Tsujii_based(eng)
print('\n')
print(Resource)