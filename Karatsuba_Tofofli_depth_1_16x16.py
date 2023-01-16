import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control


def Karatsuba_Toffoli_Depth_1(eng):
    n = 16
    a = eng.allocate_qureg(n)
    b = eng.allocate_qureg(n)
    r_a = eng.allocate_qureg(int(n / 2))
    r_b = eng.allocate_qureg(int(n / 2))
    c = eng.allocate_qureg(81)

    if (not resource_check):
        Round_constant_XOR(eng, a, 0xffee, n)
        Round_constant_XOR(eng, b, 0xeeff, n)

    new_a = []
    new_b = []
    new_r = []

    for i in range(int(n / 2)):
        CNOT | (a[i], r_a[i])
        CNOT | (a[int(n / 2) + i], r_a[i])
        CNOT | (b[i], r_b[i])
        CNOT | (b[int(n / 2) + i], r_b[i])

    new_a = karatsuba_mul_8bit(eng, a[0:8], b[0:8], c[0:27])
    new_b = karatsuba_mul_8bit(eng, a[8:16], b[8:16], c[27:54])
    new_r = karatsuba_mul_8bit(eng, r_a, r_b, c[54:81])

    # Combine
    result = []
    result = combine(eng, new_a, new_b, new_r, n)

    # Modular omitted
    if (not resource_check):
        print_state(eng, result, (2*n)-1)

def karatsuba_mul_8bit(eng, a, b, c) :

    n = 8
    r_a = eng.allocate_qureg(4)
    r_b = eng.allocate_qureg(4)

    new_a = []
    new_b = []
    new_r = []

    for i in range(int(n/2)):
        CNOT | (a[i], r_a[i])
        CNOT | (a[int(n/2)+i], r_a[i])
        CNOT | (b[i], r_b[i])
        CNOT | (b[int(n/2) + i], r_b[i])

    new_a = karatsuba_mul_4bit(eng, a[0:4], b[0:4], c[0:9])
    new_b = karatsuba_mul_4bit(eng, a[4:8], b[4:8], c[18:27])
    new_r = karatsuba_mul_4bit(eng, r_a, r_b, c[9:18])

    # Combine
    result = []
    result = combine(eng, new_a, new_b, new_r, n)

    return result

def store_operand(eng, a, b):

    c = eng.allocate_qubit()
    CNOT | (a, c)
    CNOT | (b, c)

    return c

def combine(eng, a, b, r, n):
    half_n = int(n/2)
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

def Modular_general(eng, c):

    CNOT | (c[8], c[0])
    CNOT | (c[8], c[1])
    CNOT | (c[8], c[3])
    CNOT | (c[8], c[4])

    CNOT | (c[9], c[1])
    CNOT | (c[9], c[2])
    CNOT | (c[9], c[4])
    CNOT | (c[9], c[5])

    CNOT | (c[10], c[2])
    CNOT | (c[10], c[3])
    CNOT | (c[10], c[5])
    CNOT | (c[10], c[6])

    CNOT | (c[11], c[3])
    CNOT | (c[11], c[4])
    CNOT | (c[11], c[6])
    CNOT | (c[11], c[7])

    CNOT | (c[12], c[0])
    CNOT | (c[12], c[1])
    CNOT | (c[12], c[3])
    CNOT | (c[12], c[5])
    CNOT | (c[12], c[7])

    CNOT | (c[13], c[0])
    CNOT | (c[13], c[2])
    CNOT | (c[13], c[3])
    CNOT | (c[13], c[6])

    CNOT | (c[14], c[1])
    CNOT | (c[14], c[3])
    CNOT | (c[14], c[4])
    CNOT | (c[14], c[7])


def Modular(eng, result):
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
    CNOT | (result[14], result[8])       # 8 = green
    CNOT | (result[9], result[8])

    CNOT | (result[10], result[9])            # 9 = yellow

    CNOT | (result[12], result[14]) # 14 = blue

    CNOT | (result[13], result[10])           # 10 = orange
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

def karatsuba_mul_4bit(eng, a, b, c):
    # Multiplication operands setting
    # low-part
    a0a1 = store_operand(eng, a[0], a[1])
    b0b1 = store_operand(eng, b[0], b[1])

    # middle-part
    ar0 = store_operand(eng, a[0], a[2])
    ar1 = store_operand(eng, a[1], a[3])
    br0 = store_operand(eng, b[0], b[2])
    br1 = store_operand(eng, b[1], b[3])

    # middle-inner-part
    ar0ar1 = store_operand(eng, ar0, ar1)
    br0br1 = store_operand(eng, br0, br1)

    # high-part
    a2a3 = store_operand(eng, a[2], a[3])
    b2b3 = store_operand(eng, b[2], b[3])

    # Multiplication
    # low-part
    Toffoli_gate(eng, a[0], b[0], c[0])
    Toffoli_gate(eng, a0a1, b0b1, c[1])
    Toffoli_gate(eng, a[1], b[1], c[2])

    # middle-part
    Toffoli_gate(eng, ar0, br0, c[3])
    Toffoli_gate(eng, ar0ar1, br0br1, c[4])
    Toffoli_gate(eng, ar1, br1, c[5])

    # high-part
    Toffoli_gate(eng, a[2], b[2], c[6])
    Toffoli_gate(eng, a2a3, b2b3, c[7])
    Toffoli_gate(eng, a[3], b[3], c[8])

    # Combine
    # 1

    CNOT | (c[0], c[1])
    CNOT | (c[2], c[1])

    CNOT | (c[6], c[7])
    CNOT | (c[8], c[7])

    CNOT | (c[3], c[4])
    CNOT | (c[5], c[4])

    # 2
    CNOT | (c[0], c[3])
    CNOT | (c[1], c[4])
    CNOT | (c[2], c[5])

    CNOT | (c[6], c[3])
    CNOT | (c[7], c[4])
    CNOT | (c[8], c[5])

    # 3
    CNOT | (c[2], c[3])
    CNOT | (c[6], c[5])

    output = []
    output.append(c[0])
    output.append(c[1])
    output.append(c[3])
    output.append(c[4])
    output.append(c[5])
    output.append(c[7])
    output.append(c[8])

    return output
    #print_state(eng, output, 7)

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
Karatsuba_Toffoli_Depth_1(eng)
eng.flush()
print('\n')

resource_check = 1
NCT = 0
AND_check = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
Karatsuba_Toffoli_Depth_1(eng)
print('\n')
print(Resource)