import math

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S
from projectq.backends import CircuitDrawer, ResourceCounter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control

def Karatsuba_Toffoli_Depth_1(eng) :

    n = 4
    a = eng.allocate_qureg(n) # operand a
    b = eng.allocate_qureg(n) # operand b
    c = eng.allocate_qureg(9) # result

    if (resource_check != 1):
        Round_constant_XOR(eng, a, 0xe, n)
        Round_constant_XOR(eng, b, 0xe, n)

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

    # Multiplication (Toffoli depth 1)
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

    # Combine (depth 5)
    # 1
    CNOT | (c[0], c[1])
    CNOT | (c[2], c[1])

    CNOT | (c[6], c[7])
    CNOT | (c[8], c[7])

    CNOT | (c[3], c[4])
    CNOT | (c[5], c[4])

    #2
    CNOT | (c[0], c[3])
    CNOT | (c[1], c[4])
    CNOT | (c[2], c[5])

    CNOT | (c[6], c[3])
    CNOT | (c[7], c[4])
    CNOT | (c[8], c[5])

    #3
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

    Modular(eng, output) # (x^4+x+1)

    if(resource_check != 1):
        print_state(eng, output, 4)

def Modular(eng, c):
    CNOT | (c[4], c[0]) #x^4 = x+1
    CNOT | (c[4], c[1])

    CNOT | (c[5], c[1]) #x^5 = x^2+x
    CNOT | (c[5], c[2])

    CNOT | (c[6], c[2]) #x^6 = x^3+x^2
    CNOT | (c[6], c[3])


def store_operand(eng, a, b):

    c = eng.allocate_qubit()
    CNOT | (a, c)
    CNOT | (b, c)

    return c

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
