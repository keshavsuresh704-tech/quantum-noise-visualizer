from selftest import run_selftest
from viz import run_demo, run_noise_probe


def main(state: str = "+", channel: str = "phase_flip", basis: str = "X", p: float = 0.3, shots: int = 1000, sweep: bool = True, checks: bool = False):
    if checks:
        run_selftest(backend, shots)
        run_noise_probe(backend, shots)
        return
    run_demo(backend, state, channel, basis, p, shots, sweep)