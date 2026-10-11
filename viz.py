import numpy as np
import matplotlib.pyplot as plt
from qiskit.visualization import plot_bloch_vector
from qnoise import (BASES, initial_rho, exact_channel, bloch_vector, exact_p0,
                    run_specs, backend_name)

STATE_NAMES = {"0": "|0>", "+": "|+>", "r": "|+i>", "l": "|-i>"}
CHANNEL_NAMES = {"bit_flip": "Bit flip (X error)", "phase_flip": "Phase flip (Z error)",
                 "depolarizing": "Depolarizing (random X, Y or Z)"}
BLUE, ORANGE, GREEN = "#0072B2", "#E69F00", "#009E73"      # colour-blind-safe palette
SWEEP_P = [0.0, 0.25, 0.5, 0.75, 1.0]


# ---------- input checking (the form takes free text, so be forgiving) ----------

def clean_inputs(state, channel, basis, p, shots):
    st = str(state).strip().lower().replace(" ", "")
    st = {"0": "0", "+": "+", "+i": "r", "-i": "l", "r": "r", "l": "l"}.get(st)
    if st is None:
        raise ValueError("state must be 0, + , +i or -i (you typed: " + str(state) + ")")
    ch = str(channel).strip().lower().replace(" ", "_").replace("-", "_")
    ch = {"bit_flip": "bit_flip", "bit": "bit_flip", "bitflip": "bit_flip", "x": "bit_flip",
          "phase_flip": "phase_flip", "phase": "phase_flip", "phaseflip": "phase_flip", "z": "phase_flip",
          "depolarizing": "depolarizing", "depolarising": "depolarizing", "depol": "depolarizing"}.get(ch)
    if ch is None:
        raise ValueError("channel must be bit_flip, phase_flip or depolarizing (you typed: " + str(channel) + ")")
    bs = str(basis).strip().upper()
    if bs not in BASES:
        raise ValueError("basis must be Z, X or Y (you typed: " + str(basis) + ")")
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must be between 0 and 1 (you typed: " + str(p) + ")")
    if not 10 <= shots <= 20000:
        raise ValueError("shots must be between 10 and 20000 (you typed: " + str(shots) + ")")
    return st, ch, bs, float(p), int(shots)


# ---------- exact theory ----------

def exact_after(state, channel, p):
    return exact_channel(initial_rho(state), channel, p)

def exact_curve(state, channel, basis, ps):
    rho0 = initial_rho(state)
    return [exact_p0(exact_channel(rho0, channel, float(x)), basis) for x in ps]


# ---------- sampled data from the chosen backend ----------

def sampled_p0(backend, state, channel, basis, ps, shots):
    single = "ionq" in backend_name(backend)
    if channel == "depolarizing":
        # Four fixed-Pauli circuits do not depend on p, so one run gives every p.
        specs = [(state, "pauli_" + label, 0.0, basis) for label in "IXYZ"]
        res = run_specs(backend, specs, shots=shots, single_job=single)
        values, errors = [], []
        for x in ps:
            w = [1 - x, x / 3, x / 3, x / 3]
            values.append(sum(wi * ri for wi, ri in zip(w, res)))
            var = sum(wi ** 2 * ri * (1 - ri) for wi, ri in zip(w, res)) / shots
            errors.append(float(np.sqrt(var)))
        return values, errors
    specs = [(state, channel, float(x), basis) for x in ps]
    res = run_specs(backend, specs, shots=shots, single_job=single)
    return res, [float(np.sqrt(max(r * (1 - r), 0.0) / shots)) for r in res]


# ---------- plain-language hints ----------

def blind_hint(state, channel, basis, p):
    if p == 0:
        return "p = 0 means no noise, so nothing changes. Increase p."
    rho0 = initial_rho(state)
    rho1 = exact_channel(rho0, channel, p)
    if abs(exact_p0(rho1, basis) - exact_p0(rho0, basis)) > 1e-9:
        return None
    movers = [b for b in BASES
              if abs(exact_p0(rho1, b) - exact_p0(rho0, b)) > 1e-9]
    if movers:
        return ("In basis " + basis + " this error cannot change the bars (that measurement is blind to it). "
                "Try basis " + " or ".join(movers) + ".")
    return ("This error leaves " + STATE_NAMES[state] + " completely unchanged "
            "(the state is an eigenstate of the error), so no measurement can reveal it. "
            "Try another state, for example 0 for bit flip or + for phase flip.")


# ---------- figures ----------

def show_bloch(v0, v1, title_extra):
    fig = plt.figure(figsize=(11, 5.6))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    plot_bloch_vector(list(v0), title="BEFORE noise (ideal)", ax=ax1)
    plot_bloch_vector(list(v1), title="AFTER noise (exact theory)", ax=ax2)
    fig.suptitle("Bloch sphere before and after noise: " + title_extra, fontsize=12)
    fig.text(0.27, 0.02, "<X>=%.3f  <Y>=%.3f  <Z>=%.3f" % tuple(v0), ha="center", fontsize=11)
    fig.text(0.73, 0.02, "<X>=%.3f  <Y>=%.3f  <Z>=%.3f" % tuple(v1), ha="center", fontsize=11)
    plt.show()

def show_bars(before, after, sampled, err, backend, title_extra):
    labels = ["P(0)", "P(1)"]
    groups = [
        ("Before noise (ideal)", [before, 1 - before], BLUE, ""),
        ("After noise (exact theory)", [after, 1 - after], ORANGE, "//"),
        ("After noise (measured on " + backend_name(backend) + ")", [sampled, 1 - sampled], GREEN, "xx"),
    ]
    fig, ax = plt.subplots(figsize=(9, 5.2))
    width = 0.26
    x = np.arange(2)
    for i, (name, vals, colour, hatch) in enumerate(groups):
        yerr = [err, err] if i == 2 else None
        bars = ax.bar(x + (i - 1) * width, vals, width, label=name, color=colour,
                      edgecolor="black", hatch=hatch, yerr=yerr, capsize=4)
        for rect, v in zip(bars, vals):
            ax.text(rect.get_x() + rect.get_width() / 2, min(v + 0.04, 1.1) + (0.06 if i == 2 else 0),
                    "%.3f" % v, ha="center", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylim(0, 1.2)
    ax.set_ylabel("Probability")
    ax.set_title("Measurement probabilities: " + title_extra, fontsize=11)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=1, fontsize=9)
    fig.tight_layout()
    plt.show()

def show_sweep(state, channel, basis, p, curve_ps, sampled, err, backend):
    grid = np.linspace(0, 1, 101)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(grid, exact_curve(state, channel, basis, grid), color=ORANGE, linewidth=2.5,
            label="Exact theory")
    ax.errorbar(curve_ps, sampled, yerr=err, fmt="s", color=GREEN, markeredgecolor="black",
                capsize=4, markersize=8, label="Sampled on " + backend_name(backend))
    ax.axvline(p, color="black", linestyle="--", label="Your p = %.2f" % p)
    if channel == "depolarizing":
        ax.axvline(0.75, color="gray", linestyle=":", label="p = 0.75: maximally mixed")
    ax.set_xlabel("Noise probability p")
    ax.set_ylabel("P(0) in basis " + basis)
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("Sweep over p: " + CHANNEL_NAMES[channel] + ", state " + STATE_NAMES[state]
                 + ", basis " + basis, fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


# ---------- the demo ----------

def run_demo(backend, state, channel, basis, p, shots, sweep):
    try:
        state, channel, basis, p, shots = clean_inputs(state, channel, basis, p, shots)
    except ValueError as err:
        print("Input problem:", err)
        return
    title_extra = (CHANNEL_NAMES[channel] + " on " + STATE_NAMES[state] + ", p = " + str(round(p, 3))
                   + ", measured in " + basis + " basis")
    rho0 = initial_rho(state)
    rho1 = exact_after(state, channel, p)
    v0, v1 = bloch_vector(rho0), bloch_vector(rho1)
    before, after = exact_p0(rho0, basis), exact_p0(rho1, basis)

    ps = [p]
    if sweep:
        ps += [x for x in SWEEP_P if abs(x - p) > 1e-9]
    values, errors = sampled_p0(backend, state, channel, basis, ps, shots)
    sampled, err = values[0], errors[0]

    print(title_extra)
    print("Bloch vector before: <X>=%.3f <Y>=%.3f <Z>=%.3f" % tuple(v0))
    print("Bloch vector after:  <X>=%.3f <Y>=%.3f <Z>=%.3f" % tuple(v1))
    print("P(0) before = %.3f | after (exact) = %.3f | after (sampled, %s) = %.3f +/- %.3f (1 sigma)"
          % (before, after, backend_name(backend), sampled, err))
    if abs(sampled - after) > 3 * max(err, 1.0 / shots):
        print("Note: the sampled value differs from theory by more than 3 sigma. On a noise-model "
              "simulator this is real device noise on top of our channel; on an ideal simulator it "
              "would be very unlikely.")
    hint = blind_hint(state, channel, basis, p)
    if hint:
        print("Hint:", hint)
    if channel == "depolarizing":
        print("How depolarizing is measured: four fixed-gate circuits (I, X, Y, Z) are run once and mixed "
              "with weights (1-p, p/3, p/3, p/3). The mixing is exact, so any gap from theory comes from "
              "the circuits themselves (device noise), not from shot noise in the mixing.")
    if channel == "depolarizing" and p > 0.75:
        print("Note: above p = 0.75 the Bloch vector points the opposite way. This is a valid channel: "
              "each of X, Y, Z is then more likely than doing nothing.")

    show_bloch(v0, v1, title_extra)
    show_bars(before, after, sampled, err, backend, title_extra)
    if sweep:
        order = np.argsort(ps)
        show_sweep(state, channel, basis, p, [ps[i] for i in order],
                   [values[i] for i in order], [errors[i] for i in order], backend)


# ---------- developer probe (kept from earlier tests) ----------

PROBE = [("0", "bit_flip", 0.0, "Z"), ("0", "bit_flip", 1.0, "Z"),
         ("+", "phase_flip", 0.0, "X"), ("+", "phase_flip", 1.0, "X"),
         ("r", "bit_flip", 0.0, "Y"), ("l", "bit_flip", 0.0, "Y")]

def run_noise_probe(backend, shots):
    single = "ionq" in backend_name(backend)
    sampled = run_specs(backend, PROBE, shots=shots, single_job=single)
    print("Noise probe on", backend_name(backend))
    gaps = []
    for spec, s in zip(PROBE, sampled):
        exact = exact_p0(exact_channel(initial_rho(spec[0]), spec[1], spec[2]), spec[3])
        gaps.append(abs(s - exact))
        print(spec, "exact =", round(exact, 3), "sampled =", round(s, 4), "gap =", round(abs(s - exact), 4))
    print("Mean gap:", round(sum(gaps) / len(gaps), 4), "(0.0 means noiseless)")