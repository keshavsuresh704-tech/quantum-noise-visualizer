# Quantum Noise Visualizer

**Authors:** Dylan Lee and Keshav Suresh. **Theme:** Quantum Unlocked. **License:** MIT.

This effort is supported by Qollab & IonQ.

Quantum Noise Visualizer shows what noise does to one qubit. You pick a starting state (|0>, |+>, |+i> or |-i>), a noise channel (bit flip, phase flip or depolarizing), a noise probability p, and a measurement basis (Z, X or Y). The program draws two Bloch spheres (before and after the noise), a bar chart of P(0) and P(1) comparing the ideal case, the exact theory and the values measured on a simulator, and a sweep of P(0) over p. It also tells you in plain language when your chosen measurement cannot see the error. It is built with Qiskit and runs on Qollab.

## Run it on Qollab

1. Open the project on Qollab and press **Run**.
2. Pick a simulator in the Select QPU dialog. The built-in simulator is the simplest choice. Some simulators (the IBM noise-model ones) need a free Qollab sign-in. If Run asks you to sign in, a free account is enough.
3. Fill in the form (typed values) and run:

| Field | What to type | Default |
|---|---|---|
| state | `0`, `+`, `+i` or `-i` | `+` |
| channel | `bit_flip`, `phase_flip` or `depolarizing` | `phase_flip` |
| basis | `Z`, `X` or `Y` | `X` |
| p | a number from 0 to 1 | 0.3 |
| shots | 10 to 20000 | 1000 |
| sweep | on/off: also sweep P(0) over p | on |
| checks | on/off: run the six self-checks and the noise probe instead of the demo | off |

A bad value (for example state = `7`) prints a plain message and does not crash. There is no live slider; Qollab's form takes typed values, so change p and run again.

## Run it locally

    pip install -r requirements.txt

`main.py` uses the `backend` variable that Qollab provides, so locally call the functions directly:

    python -c "from qiskit.providers.basic_provider import BasicSimulator; from viz import run_demo; run_demo(BasicSimulator(), '+', 'phase_flip', 'X', 0.3, 1000, True)"
    python -c "from qiskit.providers.basic_provider import BasicSimulator; from selftest import run_selftest; run_selftest(BasicSimulator(), 1000)"

The first command shows the demo (it opens the figures). The second runs the checks.

## Example output

Default run (phase flip on |+>, p = 0.3, measured in X, built-in simulator):

    Bloch vector before: <X>=1.000 <Y>=0.000 <Z>=0.000
    Bloch vector after:  <X>=0.400 <Y>=0.000 <Z>=0.000
    P(0) before = 1.000 | after (exact) = 0.700 | after (sampled) = 0.727 +/- 0.014 (1 sigma)

The sampled value changes a little from run to run (we saw 0.668, 0.698 and 0.727). Other runs:

| Run | Exact P(0) after | Sampled |
|---|---|---|
| Bit flip on \|0>, p = 0.3, Z basis (built-in) | 0.700 (Bloch <Z> = 0.400) | 0.699 +/- 0.015 |
| Depolarizing on \|0>, p = 0.3, Z basis (built-in) | 0.800 (Bloch <Z> = 0.600) | 0.800 +/- 0.000 |
| Bit flip on \|0>, p = 0.3, Z basis (IonQ Aria 1) | 0.700 | 0.702 +/- 0.014 |

## The physics

**Noise channels.** A noise channel changes a qubit's state at random. With probability 1 - p nothing happens. With probability p a Pauli error P is applied (X flips the bit, Z flips the phase, Y does both). In density-matrix form: rho becomes (1 - p) rho + p P rho P. The depolarizing channel applies X, Y or Z, each with probability p/3.

**What it does to the Bloch vector.** Each channel shrinks the Bloch vector along some axes:

| Channel | <X> | <Y> | <Z> |
|---|---|---|---|
| Bit flip | unchanged | x (1 - 2p) | x (1 - 2p) |
| Phase flip | x (1 - 2p) | x (1 - 2p) | unchanged |
| Depolarizing | x (1 - 4p/3) | x (1 - 4p/3) | x (1 - 4p/3) |

**Depolarizing and p = 0.75.** The factor 1 - 4p/3 is zero at p = 0.75: the state is the maximally mixed state (the centre of the sphere). Above 0.75 the factor is negative, so the Bloch vector points the opposite way, reaching -1/3 of its original vector at p = 1. This is a valid channel, not a bug: above 0.75, each of X, Y and Z is more likely than doing nothing.

**Why the measurement basis matters.** A Z measurement cannot see a phase flip on |+>: Z|+> is |->, and both give 50/50 in the Z basis. Measure in X and P(0) falls as 1 - p. The program prints a hint when your chosen basis (or state) cannot reveal the error.

**Decoherence.** Decoherence is what happens when a qubit interacts with its surroundings and loses its quantum information. The channels above are idealized, simplified models of that.

## How the noise is built: the ancilla (environment) view

For bit flip and phase flip the noise is a real circuit. A second qubit, the ancilla, acts as the environment:

1. Rotate the ancilla with RY(theta), theta = 2 arcsin(sqrt(p)), so it is |1> with probability sin^2(theta/2) = p.
2. Apply a CX (bit flip) or CZ (phase flip) from the ancilla to the data qubit.
3. Measure only the data qubit.

The joint state is cos(theta/2)|psi>|0> + sin(theta/2) P|psi>|1>. The ancilla states |0> and |1> are orthogonal, so when we ignore the ancilla there are no cross terms and the data qubit is left in exactly (1 - p) rho + p P rho P. In other words, noise is entanglement with an environment we do not look at (the Stinespring dilation).

**Depolarizing is different.** It does not use an environment qubit. Four fixed circuits (I, X, Y, Z) are run once and mixed with weights (1 - p, p/3, p/3, p/3). The mixing is exact, so the +/- value comes only from shot noise in the four circuits. It is 0.000 when each circuit gives a certain answer, as in the |0>, Z-basis example above.

The exact theory is computed with density matrices (Qiskit `quantum_info` and NumPy).

## Hardware relevance (idealized models, not claims about devices)

These are textbook models, and real devices are more complicated. Dephasing, where a qubit loses its phase information, looks like a phase-flip channel on the Bloch sphere, which is why phase errors matter for superconducting qubits. Energy relaxation (T1) is **amplitude damping**, which is not a Pauli channel (it can be approximated by one after "Pauli twirling"). Photon loss in photonic qubits is likewise not a Pauli error. So this tool shows how idealized channels behave and lets you compare them with the noise a simulated device adds on top.

## Findings (single runs on Qollab; indicative, not statistics)

The noise probe runs circuits whose ideal answer is exactly 0 or 1 and reports the gap between the sampled and ideal values. Gap of 0 means noiseless on those circuits (2000 shots each).

| Simulator | Mean gap |
|---|---|
| Built-in simulator | 0.0000 |
| IonQ Aria 1 | 0.0002 |
| IonQ Aria 2 | 0.0003 |
| IonQ Forte 1 | 0.0000 |
| IonQ Forte Enterprise 1 | 0.0000 |
| IBM Fez noise model | 0.0069 |
| IBM Kingston noise model | 0.0174 |

- The IonQ simulators' gaps are at most 0.0015 (1 to 3 counts out of 2000), so on these circuits the added noise is below what we can resolve.
- The IBM noise models add real error. In both IBM runs, the gap was larger on the circuits whose ideal answer is P(0) = 0 than on the matching circuits whose ideal answer is P(0) = 1. We saw this in one run per backend, so we do not draw conclusions from it.
- Sampling scatter on ideal and IonQ runs is consistent with shot noise (sigma is about 0.015 at p = 0.3 and 1000 shots). When a sampled value is more than 3 sigma from theory, the program says so and notes that on a noise-model simulator this is device noise on top of the injected channel.

## Verification

`selftest.py` has six checks (turn on `checks`): (1) the exact math matches the formula table, (2) the ancilla circuit equals the theory for every state, channel and basis, (3) the Y-basis rotation is right (tested on |+i> and |-i>, which have <Y> = +/-1), (4) packed circuits give the right answer for every block, (5) counts parsing (binary and hex, qubit order) and chunking are correct, (6) the four Pauli circuits mix to exactly the depolarizing channel. All six passed in Qollab. The program ran without errors in Qollab on the built-in simulator, the AWS Braket local simulator, all eight IBM noise-model simulators (Boston, Miami, Berlin, Kingston, Fez, Pittsburgh, Marrakesh, Aachen) and IonQ Aria 1.

## Engineering notes

- The built-in simulator gets about 4x slower for each extra 2-qubit block, so local runs use at most 3 blocks (6 qubits) per circuit.
- IonQ allows one job per run, so the whole sweep is packed into one circuit (12 qubits for six p values).
- The built-in simulator reports `num_qubits = None`; the code handles that. Circuits are not transpiled, because transpiling onto a 127-qubit fake backend exceeded the 24-qubit simulator limit.

## Limitations

- One qubit only. Idealized Pauli channels only; no amplitude damping.
- Depolarizing is a mixture of four circuits, not an environment model.
- No live slider (typed form values only).
- Simulators only. No real quantum hardware was used and no credits were spent.
- The IBM runs use Qollab's noise models, not real devices. The IonQ simulators showed negligible added noise on our circuits.
- Findings are single runs on small numbers of circuits.

## Future work

Amplitude-damping and coherent-error channels, two-qubit noise, repeated noise-probe runs with error bars, and a live slider.

## Files

`main.py` (form), `qnoise.py` (theory, circuits, backend runner), `selftest.py` (checks), `viz.py` (input cleaning, figures, hints, noise probe), `LICENSE`, `requirements.txt`.

## Disclosure

Built during the hackathon window with: Qiskit (with matplotlib and pylatexenc via `qiskit[visualization]`), NumPy, the Qollab Bell-state starter template, and AI assistance from Claude (Anthropic) for code, testing and documentation. The team reviewed and ran all code in Qollab. No code existed before the hackathon.

## License

MIT, see `LICENSE`. This effort is supported by Qollab & IonQ.
