# Author: Shubhodeep Chatterjee
"""Clone-level ISCN handling: tokenise a clone, resolve `idem`, relate clones.

This is deliberately not a full karyotype parser (CytoConverter, shipped with
the Mitelman download, does that). It does the one thing a copy-number
conversion flattens away: which aberrations are shared between clones of the
same investigation, which is the ordering signal.
"""
import re
from dataclasses import dataclass, field

_CELLS = re.compile(r"\[[^\]]*\]")          # [10] cell counts
_MULT = re.compile(r"x(\d+)(?:[-~]\d+)?$")    # trailing x2 / x1-4 / x2~3 multiplicity
_COUNT = re.compile(r"^[<~\d\-\s>n]+$")      # 46 | 45-46 | <2n> | 92~96
NON_EVENT = re.compile(r"^\+?(\d+[-~])?\d*(mar|dmin|r|min|hsr)\d*$|^inc$|^\?$|^idem")  # +mar, +3-4mar, 30-100dmin


def split_top(s):
    """Split on commas that are not inside parentheses."""
    out, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur))
    return [t.strip() for t in out if t.strip()]


def normalise_token(tok):
    """Return (key, multiplicity, constitutional, uncertain).

    key drops multiplicity, constitutional 'c' and question marks so that the
    same aberration written slightly differently in two clones compares equal.
    """
    tok = _CELLS.sub("", tok).strip()
    uncertain = "?" in tok
    mult = 1
    m = _MULT.search(tok)
    if m:
        mult = int(m.group(1))
        tok = tok[: m.start()]
    constitutional = tok.endswith("c") and not tok.endswith("inc")
    if constitutional:
        tok = tok[:-1]
    key = tok.replace("?", "")
    return key, mult, constitutional, uncertain


@dataclass
class Clone:
    index: int                 # 1-based order within the investigation
    text: str
    ploidy: str = ""
    uses_idem: bool = False
    idem_mult: int = 1         # K in idemxK (doubled stemline in polyploid sidelines)
    own: list = field(default_factory=list)       # aberration keys written in this clone
    abn: frozenset = frozenset()                  # full aberration set after idem expansion

    @property
    def is_normal(self):
        return not self.abn


def parse_clone(text, index=1):
    toks = split_top(_CELLS.sub("", text))
    c = Clone(index=index, text=text)
    if toks and _COUNT.match(toks[0]):
        c.ploidy = toks.pop(0)
    if toks and toks[0].startswith("idem"):
        c.uses_idem = True
        toks.pop(0)
    m = re.match(r"idemx(\d+)", c.text.split(",")[1]) if c.uses_idem and "," in c.text else None
    c.idem_mult = int(m.group(1)) if m else 1
    for t in toks:
        if re.fullmatch(r"[XY]+\??", t):          # sex-chromosome complement (may be split,
            continue                               # e.g. 47,del(X)(p22),Y,...)
        key, _, constitutional, _ = normalise_token(t)
        if constitutional or not key or NON_EVENT.match(key):
            continue
        c.own.append(key)
    return c


def parse_investigation(clone_texts):
    """Parse the clones of one investigation and expand `idem` against the stemline.

    ISCN: `idem` refers to the stemline, i.e. the first listed clone. Rules applied:
    - a sideline may drop a stemline aberration by writing it with a leading '-'
      (e.g. `-der(7)...`); honoured when the removed key exists in the stemline;
    - a numerical change opposite to one in the stemline (`+21` against `-21`)
      cancels it in a plain `idem` sideline; in a polyploid `idemxK` sideline
      (K >= 2) the doubled stemline change still dominates (two copies of -21
      plus one +21 is still a loss), so the sideline's opposite token is dropped.
    """
    clones = [parse_clone(t, i + 1) for i, t in enumerate(clone_texts)]
    if not clones:
        return clones
    stem = set(clones[0].own)
    for c in clones:
        s = set(c.own)
        if c.uses_idem:
            removed = {k for k in s if k.startswith("-") and k[1:] in stem}
            s = s - removed
            stem_c = set(stem) - {k[1:] for k in removed}
            for k in list(s):
                mnum = re.fullmatch(r"([+-])([0-9XY]+)", k)
                if not mnum:
                    continue
                opposite = ("-" if mnum.group(1) == "+" else "+") + mnum.group(2)
                if opposite in stem_c:
                    s.discard(k)
                    if c.idem_mult == 1:
                        stem_c.discard(opposite)
            s = stem_c | s
        c.abn = frozenset(s)
    return clones


def relate(clones):
    """Infer a clone forest from aberration-set containment.

    Returns dict child_index -> parent_index (None for founders) and a depth
    per clone. Parent = the abnormal clone whose set is the largest proper
    subset; ties broken by listing order. Normal clones are excluded.
    """
    ab = [c for c in clones if not c.is_normal]
    parent, depth = {}, {}
    for c in ab:
        cands = [p for p in ab if p is not c and p.abn < c.abn]
        stem = clones[0] if clones else None
        if (c.uses_idem and stem is not None and not stem.is_normal and stem is not c
                and not c.abn < stem.abn):             # never make the stemline both parent and child
            cands = cands or [stem]
        best = max(cands, key=lambda p: (len(p.abn), -p.index), default=None)
        parent[c.index] = best.index if best else None
    for c in ab:
        d, cur, seen = 0, c.index, set()
        while parent.get(cur) is not None and cur not in seen:
            seen.add(cur)
            cur = parent[cur]
            d += 1
        depth[c.index] = d
    return parent, depth


def _events(clone, key):
    ev = set()
    for a in clone.abn:
        ev |= set(key(a))
    return ev


def ordered_pairs(clones, key=lambda a: {a}, implied=False):
    """Ordering constraints implied by the clone forest from `relate`.

    For each clone with an inferred parent, every item of the parent is
    earlier than every item gained in the child. `key` maps an aberration to
    zero or more items (e.g. vocabulary events), so the same function works at
    aberration or event level. Returns a set of (earlier, later) tuples.

    implied=True also uses pairs of clones that overlap without either being
    a subset of the other (an unobserved common ancestor). Off by default:
    two clones sharing only a common late event such as -Y may well be
    independent, so these pairs are weaker evidence.
    """
    by_idx = {c.index: c for c in clones}
    parent, _ = relate(clones)
    pairs = set()
    for child, par in parent.items():
        if par is None:
            continue
        pe, ce = _events(by_idx[par], key), _events(by_idx[child], key)
        pairs |= {(e, l) for e in pe for l in ce - pe}
    if implied:
        ab = [c for c in clones if not c.is_normal]
        for a in ab:
            for b in ab:
                if a is b or a.abn <= b.abn or b.abn <= a.abn:
                    continue
                ea, eb = _events(a, key), _events(b, key)
                shared = ea & eb
                pairs |= {(e, l) for e in shared for l in eb - ea}
    return pairs
