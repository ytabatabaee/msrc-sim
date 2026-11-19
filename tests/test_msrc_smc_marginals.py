"""High-recombination marginal regression checks for the neutral SMC kernel."""
import numpy as np
import pytest
from scripts.validate_msrc_smc import _summarize_method

from msrcsim.species_tree import SpeciesTree, _Parser
from msrcsim.structured_coalescent import simulate_msc_genealogy
from msrcsim.smc import PiecewiseConstantMap, SMCConfig, simulate_msc_smc


def _balanced_species_tree(n, ne):
    branch_length = 4 * ne  # tau = 2Ne * 2, matching the high-rate case.
    def build(names, counter):
        if len(names) == 1:
            return names[0], counter
        mid = len(names) // 2
        left, counter = build(names[:mid], counter)
        right, counter = build(names[mid:], counter)
        return f"({left}:{branch_length},{right}:{branch_length})I{counter}", counter + 1

    names, _ = build([f"T{i+1}" for i in range(n)], 1)
    return SpeciesTree(names + ";", ne, 1e7)


def _root_age(newick):
    root = _Parser(newick).parse()

    def visit(node):
        if node.is_tip():
            return 0.0
        ages = [visit(child) + child.length for child in node.children]
        assert max(ages) - min(ages) < 1e-7
        node.age = float(np.mean(ages))
        return node.age

    return visit(root)


def _quartet_split(newick, chosen):
    root = _Parser(newick).parse()
    splits = []

    def visit(node):
        if node.is_tip():
            return frozenset((node.name,))
        child_sets = [visit(child) for child in node.children]
        desc = frozenset().union(*child_sets)
        side = frozenset(desc.intersection(chosen))
        other = frozenset(chosen.difference(side))
        if len(side) == len(other) == 2:
            splits.append(tuple(sorted((tuple(sorted(side)), tuple(sorted(other))))))
        return desc

    visit(root)
    assert splits, "could not resolve quartet split"
    return min(splits)


def test_quartet_summary_includes_chromosomes_with_zero_frequency():
    summaries = [{"qfreq:q": 1.0}, {}, {"qfreq:q": 0.0}]
    reference = [{"quartets": {"q": 1.0}}, {"quartets": {}}, {"quartets": {}}]
    row = next(x for x in _summarize_method(summaries, reference, "case", "method") if x["metric"] == "qfreq:q")
    assert row["method_mean"] == pytest.approx(1 / 3)
    assert row["reference_mean"] == pytest.approx(1 / 3)


@pytest.mark.parametrize("n_taxa,ne,seed", [(4, 100, 42041), (8, 10000, 82081)])
def test_high_recombination_root_age_and_quartet_marginals(n_taxa, ne, seed):
    """40 independent chromosomes vs 2,000 independent MSC trees; |z| < 3."""
    tree = _balanced_species_tree(n_taxa, ne)
    length = 10000
    rate = 20 / (4 * ne * length)
    quartet = frozenset(("T1", "T2", "T3", "T4"))
    target_split = tuple(sorted((("T1", "T2"), ("T3", "T4"))))
    linked_root = []
    linked_quartet = []
    for rep in range(40):
        result = simulate_msc_smc(
            tree,
            SMCConfig(length, PiecewiseConstantMap.constant(length, rate), seed + rep),
            np.random.default_rng(seed + rep),
        )
        widths = np.asarray([local.end - local.start for local in result.trees])
        ages = np.asarray([_root_age(local.newick) for local in result.trees])
        qcodes = np.asarray([_quartet_split(local.newick, quartet) == target_split for local in result.trees])
        linked_root.append(float(np.average(ages, weights=widths)))
        linked_quartet.append(float(np.average(qcodes, weights=widths)))

    rng = np.random.default_rng(seed + 1_000_000)
    reference_root = []
    reference_quartet = []
    for rep in range(2000):
        gene = simulate_msc_genealogy(rep, tree, rng, False)
        reference_root.append(_root_age(gene.newick))
        reference_quartet.append(_quartet_split(gene.newick, quartet) == target_split)
    reference_quartet = np.asarray(reference_quartet,dtype=bool)
    ref_qfreq = float(np.mean(reference_quartet))
    q_se = np.sqrt(np.var(linked_quartet, ddof=1) / len(linked_quartet)
                   + ref_qfreq * (1 - ref_qfreq) / len(reference_quartet))
    root_se = np.sqrt(np.var(linked_root, ddof=1) / len(linked_root)
                      + np.var(reference_root, ddof=1) / len(reference_root))
    assert abs(np.mean(linked_root) - np.mean(reference_root)) < 3 * root_se
    assert abs(np.mean(linked_quartet) - ref_qfreq) < 3 * q_se
