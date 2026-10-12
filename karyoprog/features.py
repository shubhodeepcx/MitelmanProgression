# Author: Shubhodeep Chatterjee
"""Atomic per-clone features. Vocabularies (config/vocabulary_*.json) are built on these.

Three sources:
- copy number: the union of what the ISCN text states directly (+N, -N, del(N)(q..),
  i(Nq)) and CytoConverter arm calls (karyoprog.cn) when Mitelman's syntax checker
  accepted the investigation. CytoConverter adds losses and gains hidden in
  unbalanced rearrangements (e.g. der(5)t(5;17)); the ISCN text catches interstitial
  deletions that cover less than half an arm. Features look like '+8', '-7', '5q-', '1q+'.
- sex chromosomes: from ISCN tokens only ('-Y', '-X', '+X', '+Y').
- structural and ploidy: from ISCN tokens and KCLONE modal numbers
  (e.g. 't(8;21)', 'KMT2A', 'extraPh', 'i(17q)', 'hyperdiploid51_67').
"""
import re

_T = re.compile(r"t\(([0-9XY;]+)\)\(([^)]*)\)")


def translocation_bands(tok):
    for chroms, bands in _T.findall(tok):
        cs, bs = chroms.split(";"), bands.split(";")
        if len(cs) == len(bs):
            yield dict(zip(cs, bs))


def has_t(tok, a, b, band_a=None, band_b=None):
    """True if tok contains t(a;b...) and, when bands are given in the token, the
    breakpoints start with band_a on a and band_b on b (tuples of prefixes).
    A translocation written without breakpoints is accepted."""
    if re.search(rf"t\({a};{b}[;)]", tok) is None:
        return False
    found = [bd for bd in translocation_bands(tok) if str(a) in bd and str(b) in bd]
    if not found or band_a is None:
        return True
    for bd in found:
        if a == b:                         # t(16;16), t(3;3): bands listed in order
            return True
        if bd[str(a)].startswith(band_a) and bd[str(b)].startswith(band_b):
            return True
    return False


def _band_rearranged(tok, chrom, band):
    return any(b.get(chrom, "").startswith(band) for b in translocation_bands(tok)) or \
        bool(re.match(rf"(ins|inv)\({chrom}\)\([^)]*{band}", tok))


STRUCTURAL = {
    "t(8;21)": lambda t: has_t(t, 8, 21, ("q22",), ("q22",)),
    "inv(16)": lambda t: t.startswith("inv(16)(p13") or has_t(t, 16, 16),
    "t(15;17)": lambda t: has_t(t, 15, 17, ("q2",), ("q1", "q2")),
    "t(9;22)": lambda t: has_t(t, 9, 22, ("q34",), ("q11",)),
    "extraPh": lambda t: t.startswith("+der(22)t(9;22)"),
    "t(6;9)": lambda t: has_t(t, 6, 9, ("p2",), ("q34",)),
    "inv(3)/t(3;3)": lambda t: t.startswith("inv(3)(q21") or has_t(t, 3, 3),
    "3q26": lambda t: _band_rearranged(t, "3", "q26") and not (t.startswith("inv(3)(q21") or has_t(t, 3, 3)),
    "KMT2A": lambda t: _band_rearranged(t, "11", "q23"),
    "t(12;21)": lambda t: has_t(t, 12, 21, ("p1",), ("q22",)),
    "t(1;19)": lambda t: has_t(t, 1, 19, ("q23",), ("p13",)),
    "t(8;14)": lambda t: has_t(t, 8, 14, ("q24",), ("q32",)),
    "t(11;14)": lambda t: has_t(t, 11, 14, ("q13",), ("q32",)),
    "t(4;14)": lambda t: has_t(t, 4, 14, ("p16",), ("q32",)),
    "t(14;16)": lambda t: has_t(t, 14, 16, ("q32",), ("q2",)),
    "t(14;18)": lambda t: has_t(t, 14, 18, ("q32",), ("q21",)),
    "i(17q)": lambda t: t.lstrip("+").startswith("i(17)(q10)"),
    "idic(X)(q13)": lambda t: t.startswith("idic(X)(q13"),
}

# ELN 2022: these exclude a karyotype from the "complex" category.
ELN_RECURRENT = {"t(8;21)", "inv(16)", "t(15;17)", "t(9;22)", "t(6;9)", "inv(3)/t(3;3)", "KMT2A"}


def structural_features(token):
    return {name for name, rule in STRUCTURAL.items() if rule(token)}


def iscn_cn_features(token):
    """Copy-number changes stated directly by one ISCN token."""
    m = re.fullmatch(r"([+-])(\d+)", token)
    if m:
        return {f"{m.group(1)}{m.group(2)}"}
    m = re.match(r"del\((\d+)\)\(([pq])", token)
    if m:
        return {f"{m.group(1)}{m.group(2)}-"}
    m = re.match(r"(\+?)i\((\d+)\)\(([pq])10\)", token)
    if m:
        plus, c, arm = m.groups()
        other = "q" if arm == "p" else "p"
        return {f"{c}{arm}+"} if plus else {f"{c}{arm}+", f"{c}{other}-"}
    return set()


def sex_features(token):
    return {token} if token in {"-X", "-Y", "+X", "+Y"} else set()


def ploidy_features(modal_min, modal_max):
    out = set()
    try:
        lo, hi = int(modal_min), int(modal_max)
    except (TypeError, ValueError):
        return out
    mid = (lo + hi) / 2
    if 51 <= mid <= 67:
        out.add("hyperdiploid51_67")
    if mid < 44:
        out.add("hypodiploid_lt44")
    if mid >= 81:
        out.add("near_tetraploid")
    return out


def clone_features(abn_keys, cn_calls=None, modal=(None, None), complex_threshold=3):
    """All atomic features for one clone.

    abn_keys: the clone's normalised ISCN aberration set (after idem expansion).
    cn_calls: list of CytoConverter features for this clone, or None when unavailable.
    """
    feats, struct = set(), set()
    for a in abn_keys:
        struct |= structural_features(a)
        feats |= sex_features(a)
        feats |= iscn_cn_features(a)
    if cn_calls is not None:
        feats |= set(cn_calls)
    feats |= struct
    feats |= ploidy_features(*modal)
    if len(abn_keys) >= complex_threshold and not (struct & ELN_RECURRENT):
        feats.add("complex")
    return feats


def aberration_features(token):
    """Features carried by a single aberration (used for clone-to-clone ordering)."""
    return structural_features(token) | iscn_cn_features(token) | sex_features(token)


def ploidy_range(count_token):
    """'46' -> (46, 46); '45-47' or '45~47' -> (45, 47); unparseable -> (None, None)."""
    import re as _re
    nums = [int(x) for x in _re.findall(r"\d+", count_token or "")]
    if not nums:
        return None, None
    return min(nums), max(nums)


def iscn_clone_features(clone):
    """Features of a parsed clone from its ISCN text alone (no CytoConverter)."""
    return clone_features(clone.abn, None, ploidy_range(clone.ploidy))
