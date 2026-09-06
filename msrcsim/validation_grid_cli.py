from __future__ import annotations

import argparse

from .validation_grid import default_rearrangement_fractions, run_validation_grid


def main() -> None:
    parser = argparse.ArgumentParser(description="Run v0.8.5 parameter-grid validation for the linked-spatial robustness benchmark")
    parser.add_argument("--output-dir", default="validation_grid_output")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--chrom-length", type=float, default=1_000_000.0)
    parser.add_argument("--windows", type=int, default=400)
    parser.add_argument("--block-windows", type=int, default=20)
    parser.add_argument("--kappa", type=float, default=0.25, help="Kappa used for the threshold grid")

    parser.add_argument("--tau-values", default="0.02,0.05,0.1,0.25,0.5,1.0")
    parser.add_argument("--beta-values", default="0.1,0.25,0.5,0.75,1.0")
    parser.add_argument("--rearrangement-fractions", default=",".join(f"{x:.3f}" for x in default_rearrangement_fractions()))
    parser.add_argument("--threshold-replicates", type=int, default=200)
    parser.add_argument("--threshold-bootstrap-replicates", type=int, default=500)

    parser.add_argument("--kappa-values", default="1.0,0.5,0.25,0.1,0.05")
    parser.add_argument("--kappa-rearrangement-fractions", default="0.1,0.25,0.5")
    parser.add_argument("--kappa-replicates", type=int, default=200)
    parser.add_argument("--kappa-tau", type=float, default=0.5)
    parser.add_argument("--kappa-beta", type=float, default=0.5)

    parser.add_argument("--correction-cells", default="", help="Optional semicolon-separated tau,beta,kappa cells")
    parser.add_argument("--correction-rearrangement-fractions", default="0.0,0.1,0.25,0.4,0.55,0.7")
    parser.add_argument("--correction-replicates", type=int, default=200)
    parser.add_argument("--soft-probability-mode", choices=["oracle", "noisy"], default="noisy")
    parser.add_argument("--soft-sensitivity", type=float, default=0.85)
    parser.add_argument("--soft-specificity", type=float, default=0.90)
    parser.add_argument("--soft-noise-sd", type=float, default=0.05)

    parser.add_argument("--consistency-tau", type=float, default=0.25)
    parser.add_argument("--consistency-beta", type=float, default=0.5)
    parser.add_argument("--consistency-n-blocks", default="10,25,50,100,250,500,1000")
    parser.add_argument("--consistency-replicates", type=int, default=500)
    parser.add_argument("--consistency-soft-probability-mode", choices=["oracle", "noisy"], default="oracle")
    parser.add_argument("--consistency-soft-sensitivity", type=float, default=1.0)
    parser.add_argument("--consistency-soft-specificity", type=float, default=1.0)
    parser.add_argument("--consistency-soft-noise-sd", type=float, default=0.0)

    args = parser.parse_args()
    out = run_validation_grid(args)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
