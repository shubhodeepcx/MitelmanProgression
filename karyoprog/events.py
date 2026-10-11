"""Provisional event vocabulary: map one normalised ISCN aberration to events.

PROVISIONAL (week 1). The frozen vocabulary is a weeks 2-3 deliverable; see
docs/event_vocabulary.md. Copy-number events here are read from the ISCN token
only; losses through unbalanced translocations (e.g. der(5)t(5;17)) are not yet
counted and will come from the CytoConverter segments.
"""
import re

_T = re.compile(r"t\(([0-9XY;]+)\)\(([^)]*)\)")


def _translocation_bands(tok):
    """Yield {chrom: band} for every t(...)(...) inside a token."""
    for chroms, bands in _T.findall(tok):
        cs, bs = chroms.split(";"), bands.split(";")
        if len(cs) == len(bs):
            yield dict(zip(cs, bs))


def _has_t(tok, a, b):
    return re.search(rf"t\({a};{b}[;)]", tok) is not None


def _num(tok, sign, chrom):
    return tok == f"{sign}{chrom}"


RULES = {
    # recurrent balanced rearrangements (one event each)
    "t(8;21)": lambda t: _has_t(t, 8, 21),
    "inv(16)": lambda t: t.startswith("inv(16)") or _has_t(t, 16, 16),
    "t(15;17)": lambda t: _has_t(t, 15, 17),
    "t(9;22)": lambda t: _has_t(t, 9, 22),
    "t(6;9)": lambda t: _has_t(t, 6, 9),
    "inv(3)/t(3;3)": lambda t: t.startswith("inv(3)(q21") or _has_t(t, 3, 3),
    "KMT2A(11q23)": lambda t: any(b.get("11", "").startswith("q23") for b in _translocation_bands(t))
                              or bool(re.match(r"(ins|inv)\(11\)\([^)]*q23", t)),
    # whole-chromosome and arm-level copy number
    "-5/del(5q)": lambda t: _num(t, "-", 5) or t.startswith("del(5)(q"),
    "-7": lambda t: _num(t, "-", 7),
    "del(7q)": lambda t: t.startswith("del(7)(q"),
    "+8": lambda t: _num(t, "+", 8),
    "-17/17p-": lambda t: _num(t, "-", 17) or t.startswith("del(17)(p") or t.startswith("i(17)(q10)"),
    "del(20q)": lambda t: t.startswith("del(20)(q"),
    "del(9q)": lambda t: t.startswith("del(9)(q"),
    "+21": lambda t: _num(t, "+", 21),
    "+22": lambda t: _num(t, "+", 22),
    "+13": lambda t: _num(t, "+", 13),
    "+11": lambda t: _num(t, "+", 11),
    "+4": lambda t: _num(t, "+", 4),
    "-Y": lambda t: _num(t, "-", "Y"),
    "-X": lambda t: _num(t, "-", "X"),
}

# Balanced rearrangements that, per ELN 2017/2022, exclude a karyotype from
# the "complex" category.
ELN_RECURRENT = {"t(8;21)", "inv(16)", "t(15;17)", "t(9;22)", "t(6;9)",
                 "inv(3)/t(3;3)", "KMT2A(11q23)"}


def events_of(token):
    return [name for name, rule in RULES.items() if rule(token)]


def is_complex(abn_set, threshold=3):
    """ELN-style complex karyotype: >=3 aberrations, no recurrent balanced event."""
    evs = {e for a in abn_set for e in events_of(a)}
    return len(abn_set) >= threshold and not (evs & ELN_RECURRENT)
