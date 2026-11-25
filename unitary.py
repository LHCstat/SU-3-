"""U(3) and SU(3) via complex plane rotations. All angles in radians."""
import argparse
import numpy as np

def rotation(i, j, theta, phi=0):
    """Zero-based coordinate indices; exp(theta G_ij)."""
    if not 0 <= i < j < 3:
        raise ValueError("Require 0 <= i < j < 3")
    G = np.zeros((3, 3), dtype=complex)
    G[i, j] = np.exp(-1j * phi)
    G[j, i] = -np.exp(1j * phi)
    return np.eye(3) + np.sin(theta)*G + (1-np.cos(theta))*(G@G)

def unitary3(theta12=0, theta13=0, theta23=0,
             phi12=0, phi13=0, phi23=0,
             alpha1=0, alpha2=0, alpha3=0, *, special=False):
    """Generate U(3); with special=True generate SU(3), overriding alpha3."""
    if special:
        alpha3 = -alpha1-alpha2
    D = np.diag(np.exp(1j*np.array([alpha1, alpha2, alpha3])))
    return (D @ rotation(1, 2, theta23, phi23)
              @ rotation(0, 2, theta13, phi13)
              @ rotation(0, 1, theta12, phi12))

def check(U):
    return float(np.max(np.abs(U.conj().T@U-np.eye(3)))), np.linalg.det(U)

def self_test():
    rng = np.random.default_rng(42)
    worst = 0.
    for special in (False, True):
        for _ in range(1000):
            p = rng.uniform(0, 2*np.pi, 9)
            U = unitary3(*p, special=special)
            err, det = check(U)
            expected = 1 if special else np.exp(1j*sum(p[6:]))
            assert err < 1e-12 and abs(det-expected) < 1e-12
            worst = max(worst, err)
    print(f"2000 matrices passed; maximum error = {worst:.3e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--angles", nargs=9, type=float, help="theta12 theta13 theta23 phi12 phi13 phi23 alpha1 alpha2 alpha3")
    parser.add_argument("--degrees", action="store_true")
    parser.add_argument("--su3", action="store_true")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.test:
        self_test()
    else:
        if args.angles is None:
            rng = np.random.default_rng(args.seed)
            p = np.r_[rng.uniform(0,np.pi/2,3),rng.uniform(0,2*np.pi,6)]
        else:
            p = np.deg2rad(args.angles) if args.degrees else args.angles
        U = unitary3(*p, special=args.su3)
        err, det = check(U)
        np.set_printoptions(precision=8, suppress=True)
        print("U =", U, "\nU†U =", U.conj().T@U, "\ndet(U) =", det, "\nmax error =", err)
