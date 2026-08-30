from __future__ import annotations

from pathlib import Path
from typing import Any
import argparse

import yaml

from .analytic import TOPOLOGY_NAMES
from .hybridization_match import fit_hybridization_to_q
from .spatial_compare_io import load_spatial_model_run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fit pulse-hybridization parameters to match an MSRC marginal quartet vector")
    parser.add_argument("--msrc-dir", required=True, help="Directory containing MSRC spatial outputs")
    parser.add_argument("--major-topology", required=True, help="Major parental topology: 12|34, 13|24, 14|23, or 0/1/2")
    parser.add_argument("--introgressed-topology", required=True, help="Introgressed parental topology: 12|34, 13|24, 14|23, or 0/1/2")
    parser.add_argument("--output", required=True, help="Output YAML path")
    parser.add_argument("--max-branch-length", type=float, default=20.0, help="Upper bound for fitted MSC branch lengths")
    parser.add_argument("--chromosome-length-bp", type=float, help="Override chromosome length recorded in the generated config block")
    parser.add_argument("--num-loci", type=int, default=5000, help="Default loci count recorded in the generated config block")
    parser.add_argument("--h", type=float, default=200.0, help="Default hybridization age recorded in the generated config block")
    parser.add_argument("--r", type=float, default=1e-8, help="Default recombination rate recorded in the generated config block")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        msrc_run = load_spatial_model_run(args.msrc_dir, "msrc")
        match = fit_hybridization_to_q(
            msrc_run.marginal_q,
            args.major_topology,
            args.introgressed_topology,
            max_branch_length=args.max_branch_length,
        )
        length = float(args.chromosome_length_bp or msrc_run.chromosome_length_bp)
        expected_intro_len = None if args.h <= 0 or args.r <= 0 or match.gamma >= 1.0 else 1.0 / (args.h * args.r * max(1e-15, 1.0 - match.gamma))
        data: dict[str, Any] = {
            "mode": "matched_pulse_hybridization",
            "msrc_dir": str(msrc_run.run_dir),
            "target_q": list(match.target_q),
            "fitted_q": list(match.fitted_q),
            "l1_error": match.l1_error,
            "l2_error": match.l2_error,
            "fit_success": match.success,
            "fit_message": match.message,
            "hybridization": {
                "gamma": match.gamma,
                "t_major": match.t_major,
                "t_introgressed": match.t_introgressed,
                "major_topology": match.major_topology,
                "major_topology_label": TOPOLOGY_NAMES[match.major_topology],
                "introgressed_topology": match.introgressed_topology,
                "introgressed_topology_label": TOPOLOGY_NAMES[match.introgressed_topology],
                "generations_since_hybridization": float(args.h),
                "recombination_rate": float(args.r),
                "expected_mean_introgressed_tract_length_bp": expected_intro_len,
            },
            "simulation_template": {
                "chromosome": {
                    "length_bp": length,
                    "num_loci": int(args.num_loci),
                    "locus_positions": {"mode": "evenly_spaced"},
                },
                "windows": {"loci_per_window": 50, "step_loci": 10},
                "spatial_statistics": {"lags_bp": [10000, 100000, 1000000, 5000000]},
            },
        }
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(yaml.safe_dump(data, sort_keys=False))
        print(f"Fitted gamma={match.gamma:.6g}, t_major={match.t_major:.6g}, t_intro={match.t_introgressed:.6g}")
        print(f"Target q={match.target_q}; fitted q={match.fitted_q}; L1={match.l1_error:.6g}; L2={match.l2_error:.6g}")
        print(f"Wrote {output}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
