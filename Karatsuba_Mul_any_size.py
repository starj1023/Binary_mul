import math
from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control

def Karatsuba_Toffoli_Depth_1(eng) :

    n = 127
    a = eng.allocate_qureg(n)
    b = eng.allocate_qureg(n)

    if (resource_check != 1):
        Round_constant_XOR(eng, a, 0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff, n)
        Round_constant_XOR(eng, b, 0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff, n)

    result = []
    result= recursive_karatsuba(eng, a, b, n)

    if (resource_check != 1):
        print('\nProduct: ', end='')
        print_state(eng, result, 2 * n - 1)

def recursive_karatsuba(eng, a, b, n): #n=4

    if(n==1):
        c = eng.allocate_qubit()
        Toffoli_gate(eng, a, b, c)

        return c

    c_len = 3**math.log(n, 2) #9 #3
    r_low = n//2    #2 #1

    if(n%2!=0):
        r_low = r_low +1 # n=3 -> 2, n=4 -> 2

    r_a = []
    r_b = []

    # Provide rooms and prepare operands
    r_a = room(eng, r_low) #2qubits for r
    r_b = room(eng, r_low) #2qubits for r

    # r_a = a_low + a_high
    copy(eng, a[0:r_low], r_a, r_low)
    copy(eng, a[r_low:n], r_a, n//2)

    # r_b = b_low + b_high
    copy(eng, b[0:r_low], r_b, r_low)
    copy(eng, b[r_low:n], r_b, n//2)

    # upper-part setting
    if(r_low == 1):
        c = eng.allocate_qureg(3)
        Toffoli_gate(eng, a[0], b[0], c[0])
        Toffoli_gate(eng, a[1], b[1], c[2])
        CNOT | (c[0], c[1])
        CNOT | (c[2], c[1])
        Toffoli_gate(eng, r_a, r_b, c[1])

        return c

    c_a = []
    c_b = []
    c_r = []

    c_a = recursive_karatsuba(eng, a[0:r_low], b[0:r_low], r_low)# 2 qubits     # 0~2
    c_b = recursive_karatsuba(eng, a[r_low:n], b[r_low:n], n//2)#2 qubits        # 3~5
    c_r = recursive_karatsuba(eng, r_a[0:r_low], r_b[0:r_low], r_low) #2qubits  # 6~8

    result = []
    result = combine(eng, c_a, c_b, c_r, n)

    return result

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

    ###### 7-bit
    for i in range(7):
        CNOT | (a[i], r[i])
    for i in range(5):
        CNOT | (b[i], r[i])

    for i in range(3):
        CNOT | (a[4+i], r[i])
    for i in range(3):
        CNOT | (b[i], r[4+i])

    out = []
    for i in range(4):
        out.append(a[i])
    for i in range(7):
        out.append(r[i])
    out.append(b[3])
    out.append(b[4])

    return out

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

resource_check = 1
NCT = 0
AND_check = 0
Resource = ResourceCounter()
eng = MainEngine(Resource)
Karatsuba_Toffoli_Depth_1(eng)
print('\n')
print(Resource)