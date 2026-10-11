import numpy as np
from qiskit.quantum_info import DensityMatrix, partial_trace
from qnoise import (STATES, CHANNELS, BASES, initial_rho, exact_channel, bloch_vector,
                    exact_p0, pauli_conjugate, build_circuit, pack_circuit, marginal_p0, plan_chunks,
                    run_specs, backend_name)

def exact_value(state, channel, p, basis):
    return exact_p0(exact_channel(initial_rho(state), channel, p), basis)

def circuit_p0(state, channel, p, basis):
    qc = build_circuit(state, channel, p, basis, measure=False)
    reduced = partial_trace(DensityMatrix(qc), [1])
    return float(np.real(reduced.probabilities()[0]))

def check_exact_math():
    for p in (0.0, 0.3, 0.75, 1.0):
        for state in STATES:
            v0 = bloch_vector(initial_rho(state))
            factors = {"bit_flip": np.array([1, 1 - 2 * p, 1 - 2 * p]),
                       "phase_flip": np.array([1 - 2 * p, 1 - 2 * p, 1]),
                       "depolarizing": np.full(3, 1 - 4 * p / 3)}
            for ch, f in factors.items():
                got = bloch_vector(exact_channel(initial_rho(state), ch, p))
                assert np.allclose(got, v0 * f), (state, ch, p)
    print("Check 1 passed: exact math matches the formula table")

def check_circuit_equals_theory():
    for state in STATES:
        for ch in CHANNELS:
            for basis in BASES:
                for p in (0.0, 0.1, 0.3, 0.5, 0.9, 1.0):
                    assert np.isclose(circuit_p0(state, ch, p, basis),
                                      exact_value(state, ch, p, basis)), (state, ch, basis, p)
    print("Check 2 passed: ancilla circuit equals theory (all states, all bases)")

def check_y_basis():
    assert np.isclose(circuit_p0("r", "bit_flip", 0.0, "Y"), 1.0)
    assert np.isclose(circuit_p0("l", "bit_flip", 0.0, "Y"), 0.0)
    assert np.isclose(circuit_p0("r", "bit_flip", 0.3, "Y"), 0.7)
    print("Check 3 passed: Y-basis rotation is correct (tested on states with <Y> = +/-1)")

def check_pauli_variants():
    for state in STATES:
        for basis in BASES:
            ps = []
            for label in "IXYZ":
                qc = build_circuit(state, "pauli_" + label, 0.0, basis, measure=False)
                reduced = partial_trace(DensityMatrix(qc), [1])
                got = float(np.real(reduced.probabilities()[0]))
                want = exact_p0(pauli_conjugate(initial_rho(state), label), basis)
                assert np.isclose(got, want), (state, label, basis)
                ps.append(got)
            for p in (0.0, 0.3, 0.75, 1.0):
                mix = (1 - p) * ps[0] + (p / 3) * (ps[1] + ps[2] + ps[3])
                assert np.isclose(mix, exact_value(state, "depolarizing", p, basis)), (state, basis, p)
    print("Check 6 passed: the four Pauli circuits mix to exactly the depolarizing channel")

def check_packer():
    specs = [("0", "bit_flip", 0.3, "Z"), ("+", "phase_flip", 0.6, "X"),
             ("r", "bit_flip", 0.2, "Y"), ("+", "phase_flip", 0.9, "Z")]
    n = 2 * len(specs)
    dm = DensityMatrix(pack_circuit(specs, measure=False))
    for i, (state, ch, p, basis) in enumerate(specs):
        others = [j for j in range(n) if j != 2 * i]
        got = float(np.real(partial_trace(dm, others).probabilities()[0]))
        assert np.isclose(got, exact_value(state, ch, p, basis)), (i, specs[i])
    print("Check 4 passed: packed circuit gives the right answer for every block")

def check_marginals_and_chunks():
    demo = {"110": 30, "001": 70}
    assert np.isclose(marginal_p0(demo, 0), 0.3)
    assert np.isclose(marginal_p0(demo, 2), 0.7)
    assert np.isclose(marginal_p0({"0x6": 30, "0x1": 70}, 0, 3), 0.3)
    assert [len(c) for c in plan_chunks(list(range(7)), 5)] == [2, 2, 2, 1]
    assert [len(c) for c in plan_chunks(list(range(7)), None)] == [3, 3, 1]
    assert [len(c) for c in plan_chunks(list(range(7)), 127)] == [3, 3, 1]
    try:
        plan_chunks([1], 1)
        raise AssertionError("should have rejected a 1-qubit backend")
    except ValueError:
        pass
    print("Check 5 passed: marginals (binary and hex) and chunking correct")

def run_selftest(backend=None, shots=2000):
    check_exact_math()
    check_circuit_equals_theory()
    check_y_basis()
    check_packer()
    check_marginals_and_chunks()
    check_pauli_variants()
    if backend is None:
        print("No backend smoke test requested.")
        return
    specs = [("0", "bit_flip", 0.3, "Z"), ("+", "phase_flip", 0.3, "X"),
             ("+", "phase_flip", 0.3, "Z")]
    single = "ionq" in backend_name(backend)
    sampled = run_specs(backend, specs, shots=shots, single_job=single)
    print("Backend:", backend_name(backend), "| single job:", single)
    for spec, s in zip(specs, sampled):
        print(spec, "exact =", round(exact_value(*spec), 3), "sampled =", round(s, 3))