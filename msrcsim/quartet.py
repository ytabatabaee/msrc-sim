from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Iterable


@dataclass
class ParsedNewickNode:
    name: str = ""
    length: float = 0.0
    children: list["ParsedNewickNode"] = field(default_factory=list)

    def is_tip(self) -> bool:
        return not self.children


class _Parser:
    def __init__(self, text: str):
        self.text = text.strip().rstrip(";")
        self.i = 0

    def parse(self) -> ParsedNewickNode:
        node = self.node()
        if self.i != len(self.text):
            raise ValueError(f"Unexpected Newick text at {self.text[self.i:]}")
        return node

    def node(self) -> ParsedNewickNode:
        if self.text[self.i] == "(":
            self.i += 1
            children = [self.node()]
            while self.i < len(self.text) and self.text[self.i] == ",":
                self.i += 1
                children.append(self.node())
            if self.i >= len(self.text) or self.text[self.i] != ")":
                raise ValueError("Malformed Newick")
            self.i += 1
            name = self.label()
            length = self.length()
            return ParsedNewickNode(name, length, children)
        name = self.label()
        if not name:
            raise ValueError("Tip names are required")
        return ParsedNewickNode(name, self.length(), [])

    def label(self) -> str:
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in ":,()":
            self.i += 1
        return self.text[start:self.i].strip()

    def length(self) -> float:
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            start = self.i
            while self.i < len(self.text) and self.text[self.i] not in ",()":
                self.i += 1
            value = float(self.text[start:self.i])
            if value < -1e-12:
                raise ValueError("Newick branch lengths must be nonnegative")
            return value
        return 0.0


def parse_newick(text: str) -> ParsedNewickNode:
    return _Parser(text).parse()


def leaf_names(node: ParsedNewickNode) -> list[str]:
    if node.is_tip():
        return [node.name]
    out: list[str] = []
    for child in node.children:
        out.extend(leaf_names(child))
    return out


def validate_newick_taxa(text: str, expected_taxa: Iterable[str]) -> None:
    root = parse_newick(text)
    observed = leaf_names(root)
    expected = list(expected_taxa)
    if sorted(observed) != sorted(expected):
        raise ValueError(f"Expected taxa {sorted(expected)}, observed {sorted(observed)}")
    if len(observed) != len(set(observed)):
        raise ValueError("Newick contains duplicate taxon labels")


def _prune(node: ParsedNewickNode, keep: set[str]) -> ParsedNewickNode | None:
    if node.is_tip():
        return ParsedNewickNode(node.name, node.length, []) if node.name in keep else None
    children = [child for child in (_prune(c, keep) for c in node.children) if child is not None]
    if not children:
        return None
    if len(children) == 1:
        child = children[0]
        child.length += node.length
        return child
    return ParsedNewickNode(node.name, node.length, children)


def _to_newick(node: ParsedNewickNode, is_root: bool = True) -> str:
    if node.is_tip():
        text = node.name
    else:
        text = "(" + ",".join(_to_newick(child, False) for child in node.children) + ")" + node.name
    if not is_root:
        text += f":{max(0.0, node.length):.8f}"
    return text


def prune_newick_to_taxa(newick: str, taxa: Iterable[str]) -> str:
    keep = set(taxa)
    if len(keep) != 4:
        raise ValueError("Quartet pruning requires exactly four distinct taxa")
    pruned = _prune(parse_newick(newick), keep)
    if pruned is None:
        raise ValueError("None of the requested taxa occur in the Newick tree")
    observed = set(leaf_names(pruned))
    if observed != keep:
        raise ValueError(f"Missing quartet taxa: {sorted(keep - observed)}")
    return _to_newick(pruned) + ";"


def _descendant_sets(node: ParsedNewickNode) -> list[frozenset[str]]:
    if node.is_tip():
        return [frozenset([node.name])]
    out = []
    leaves = frozenset(leaf_names(node))
    out.append(leaves)
    for child in node.children:
        out.extend(_descendant_sets(child))
    return out


def classify_quartet_newick(newick: str, taxa: Iterable[str]) -> int:
    ordered = list(taxa)
    if len(ordered) != 4 or len(set(ordered)) != 4:
        raise ValueError("Quartet classification requires exactly four distinct taxa")
    pruned = parse_newick(prune_newick_to_taxa(newick, ordered))
    splits = {s for s in _descendant_sets(pruned) if len(s) == 2}
    pairs = [
        frozenset((ordered[0], ordered[1])),
        frozenset((ordered[0], ordered[2])),
        frozenset((ordered[0], ordered[3])),
    ]
    full = frozenset(ordered)
    for idx, pair in enumerate(pairs):
        if pair in splits or (full - pair) in splits:
            return idx
    # Star-like or unresolved trees are assigned to topology 0 for compatibility
    # with the legacy fallback classifier.
    return 0


def summarize_quartet(gene_trees: Iterable[str], taxa: Iterable[str]) -> dict[str, object]:
    ordered = list(taxa)
    counts = [0, 0, 0]
    for tree in gene_trees:
        counts[classify_quartet_newick(tree, ordered)] += 1
    total = sum(counts)
    q = [count / total if total else float("nan") for count in counts]
    return {
        "taxa": ordered,
        "n1": counts[0],
        "n2": counts[1],
        "n3": counts[2],
        "q1": q[0],
        "q2": q[1],
        "q3": q[2],
    }


def all_quartets(taxa: Iterable[str]) -> list[tuple[str, str, str, str]]:
    return [tuple(c) for c in combinations(list(taxa), 4)]
