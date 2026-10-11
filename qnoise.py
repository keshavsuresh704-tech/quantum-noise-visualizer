import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import DensityMatrix, Operator, Pauli

STATES = ("0", "+", "r", "l")
CHANNELS = ("bit_flip", "phase_flip")
BASES = ("Z", "X", "Y")
PAULI_VARIANTS = ("pauli_I", "pauli_X", "pauli_Y", "pauli_Z")   # fixed gates, used to sample depolarizing

def initial_rho(state):
    return DensityMatrix.from_label(state)

def pauli_conjugate(rho, label):
    return rho.evolve(Operator(Pauli(label)))

def exact_channel(rho, channel, p):
    paulis = {"bit_flip": ["X"], "phase_flip": ["Z"],
              "depolarizing": ["X", "Y", "Z"]}[channel]
    weight = p / len(paulis)
    data = (1 - p) * rho.data
    for label in paulis:
        data = data + weight * pauli_conjugate(rho, label).data
    return DensityMatrix(data)

def bloch_vector(rho):
    return np.real([rho.expectation_value(Pauli(s)) for s in "XYZ"])

def exact_p0(rho, basis):
    return (1 + bloch_vector(rho)["XYZ".index(basis)]) / 2
  
def add_block(qc, data, anc, clbit, state, channel, p, basis, measure=True):
    if state not in STATES:
        raise ValueError("state must be one of " + str(STATES))
    if channel not in CHANNELS and channel not in PAULI_VARIANTS:
        raise ValueError("channel must be 'bit_flip' or 'phase_flip' (or a pauli_ variant)")
    if basis not in BASES:
        raise ValueError("basis must be 'Z', 'X' or 'Y'")
    if not 0 <= p <= 1:
        raise ValueError("p must be between 0 and 1")

    if state in ("+", "r", "l"):
        qc.h(data)
    if state == "r":
        qc.s(data)
    elif state == "l":
        qc.sdg(data)
    if channel in PAULI_VARIANTS:
        gate = channel[-1]               # a fixed Pauli, no environment qubit needed
        if gate == "X":
            qc.x(data)
        elif gate == "Y":
            qc.y(data)
        elif gate == "Z":
            qc.z(data)
    else:
        theta = 2 * np.arcsin(np.sqrt(np.clip(p, 0.0, 1.0)))
        qc.ry(theta, anc)
        if channel == "bit_flip":
            qc.cx(anc, data)
        else:
            qc.cz(anc, data)
    qc.barrier(data, anc)
    if basis == "X":
        qc.h(data)
    elif basis == "Y":
        qc.sdg(data)
        qc.h(data)
    if measure:
        qc.measure(data, clbit)

def build_circuit(state, channel, p, basis, measure=True):
    qc = QuantumCircuit(2, 1)
    add_block(qc, 0, 1, 0, state, channel, p, basis, measure)
    return qc

def pack_circuit(specs, measure=True):
    k = len(specs)
    if k == 0:
        raise ValueError("need at least one spec")
    qc = QuantumCircuit(2 * k, k)
    for i, (state, channel, p, basis) in enumerate(specs):
        add_block(qc, 2 * i, 2 * i + 1, i, state, channel, p, basis, measure)
    return qc

def marginal_p0(counts, clbit, nbits=None):
    total = 0
    zeros = 0
    for key, n in counts.items():
        text = str(key).replace(" ", "")
        if text.startswith("0x"):
            if nbits is None:
                raise ValueError("hex counts keys need nbits")
            text = format(int(text, 16), "0" + str(nbits) + "b")
        if nbits is not None:
            text = text.zfill(nbits)
        total += n
        if text[::-1][clbit] == "0":
            zeros += n
    return zeros / total

LOCAL_MAX_BLOCKS = 3   # BasicSimulator slows ~4x per extra block; keep local circuits <= 6 qubits

def backend_qubits(backend):
    n = backend.num_qubits
    return n() if callable(n) else n          # may be None (built-in simulator)

def backend_name(backend):
    n = getattr(backend, "name", "")
    n = n() if callable(n) else n
    return str(n).lower()

def plan_chunks(specs, num_qubits, max_blocks=LOCAL_MAX_BLOCKS):
    if num_qubits is not None and num_qubits < 2:
        raise ValueError("backend has " + str(num_qubits) + " qubit(s); need at least 2")
    per_circuit = max_blocks if num_qubits is None else min(max_blocks, num_qubits // 2)
    specs = list(specs)
    return [specs[i:i + per_circuit] for i in range(0, len(specs), per_circuit)]

def run_specs(backend, specs, shots=1000, single_job=False):
    num_qubits = backend_qubits(backend)
    if single_job:
        if num_qubits is not None and 2 * len(specs) > num_qubits:
            raise ValueError("need " + str(2 * len(specs)) + " qubits in one job, backend has "
                             + str(num_qubits) + "; use fewer p values")
        chunks = [list(specs)]
    else:
        chunks = plan_chunks(specs, num_qubits)
    results = []
    for chunk in chunks:
        qc = pack_circuit(chunk)
        counts = backend.run(qc, shots=shots).result().get_counts()
        for i in range(len(chunk)):
            results.append(marginal_p0(counts, i, len(chunk)))
    return results