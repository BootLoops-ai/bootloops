"""emitall.expr — the claims-spec field-expression language.

A field expression is a short Python-syntax string evaluated by a whitelisted
AST walker (never `eval`). It can reference receipt values only through the
functions below — there is no attribute access, no subscripting, no name
binding, no lambdas — so a spec cannot execute anything beyond arithmetic and
the registered accessors.

Accessors
---------
    f(alias, pointer)          value at a JSON pointer inside receipt `alias`
                               (segments split on "/", numeric segment = list
                               index, "~0"/"~1" escapes per RFC 6901)
    vals(alias, pattern)       list of values; pattern segments may be "*"
                               (every key/index, document order) or "{a|b|c}"
                               (alternation, listed order)
    rows(alias, pointer="")    list-of-rows at pointer ("" = whole receipt;
                               a CSV receipt is already a list of dict rows)
    cw(alias, pointer, conds)  count-where: rows matching ALL conditions.
                               conds is a JSON string: a list of
                               [field, op, value] or [field, op, value, coerce]
                               with op in ==, !=, >, >=, <, <=, in, not_in;
                               coerce in "int"/"float" applies to the ROW field
                               before comparing; value may be
                               {"$ref": "alias:/pointer"} to compare against
                               (or test membership in) another receipt's value.
                               Conditions are evaluated in listed order with
                               short-circuit AND, so a coercion can be guarded
                               by a preceding not_in/!= condition.
    joinrows(alias, pointer, template, sep)
                               sep.join(template.format(**row) for each row)

Helpers (pure): len sum min max abs int float str round bool,
    where(xs, op, v)  -> filtered list        split(s, sep, i) -> str piece
    join(xs, sep)     -> sep-joined strings   jsondump(x) -> json.dumps(
                                                             x, sort_keys=True)
    zipjoin(xs, ys, template, sep)            row-wise zip of two list fields:
                               sep.join(template.format(a=x, b=y)); LOUD on
                               length mismatch
    frac(x)           -> exact Fraction from "3/7", "0.1", "1e-6", int, or a
                               float's shortest decimal — Fractions flow
                               through the arithmetic operators, str() prints
                               "p/q", float() converts; for exact interval
                               widths

Condition/where ops: ==, !=, >, >=, <, <=, in, not_in, and
contains / not_contains — substring (or element) predicate testing
`value in element`, e.g. where(L, 'contains', 'FY20') or a cw condition
["field", "contains", "restated"]; LOUD on a non-container element.

Operators: + - * / // % **, unary -, comparisons (==, !=, <, <=, >, >=,
in, not in), `and`/`or`/`not`, list/tuple/dict literals. Comparisons yield
True/False and may be summed (Python bool arithmetic) — that is how
"count of groups with a unique optimum" style checks are written.

Errors are LOUD: unknown alias/name/function, missing pointer segment, or
disallowed syntax raises ExprError (ReceiptMissing for an alias whose file
was absent at load). The engine turns these into findings, never silence.
No bare `assert` anywhere (survives `python -O`).
"""
import ast
import json
from fractions import Fraction


class ExprError(ValueError):
    """A field expression could not be evaluated against the receipts."""


class ReceiptMissing(ExprError):
    """The expression touched a receipt whose file was absent at load time."""

    def __init__(self, alias):
        super().__init__(f"receipt {alias!r} missing on disk")
        self.alias = alias


_MISSING = object()  # sentinel stored by the engine for absent receipt files


class Env:
    """alias -> loaded receipt object (or the missing sentinel)."""

    MISSING = _MISSING

    def __init__(self, receipts):
        self.receipts = receipts

    def get(self, alias):
        if alias not in self.receipts:
            raise ExprError(f"unknown receipt alias {alias!r}")
        obj = self.receipts[alias]
        if obj is _MISSING:
            raise ReceiptMissing(alias)
        return obj


# ---------------------------------------------------------------- pointers

def _unescape(seg):
    return seg.replace("~1", "/").replace("~0", "~")


def _child(obj, seg):
    if isinstance(obj, list):
        try:
            return obj[int(seg)]
        except (ValueError, IndexError) as e:
            raise ExprError(f"bad list index {seg!r}: {e}") from None
    if isinstance(obj, dict):
        if seg not in obj:
            raise ExprError(f"missing key {seg!r} (have {sorted(obj)[:8]}...)")
        return obj[seg]
    raise ExprError(f"cannot descend into {type(obj).__name__} with {seg!r}")


def resolve(obj, pointer):
    """Value at a JSON pointer ("" or "/" = the whole object)."""
    if not isinstance(pointer, str):
        raise ExprError(f"pointer must be a string, got {pointer!r}")
    if pointer in ("", "/"):
        return obj
    for seg in pointer.lstrip("/").split("/"):
        obj = _child(obj, _unescape(seg))
    return obj


def resolve_pattern(obj, pattern):
    """Flat list of values matched by a pointer pattern with * / {a|b|c}."""
    if pattern in ("", "/"):
        return [obj]
    out = [obj]
    for seg in pattern.lstrip("/").split("/"):
        nxt = []
        for o in out:
            if seg == "*":
                if isinstance(o, dict):
                    nxt.extend(o.values())
                elif isinstance(o, list):
                    nxt.extend(o)
                else:
                    raise ExprError(f"'*' cannot expand {type(o).__name__}")
            elif seg.startswith("{") and seg.endswith("}"):
                for k in seg[1:-1].split("|"):
                    nxt.append(_child(o, _unescape(k)))
            else:
                nxt.append(_child(o, _unescape(seg)))
        out = nxt
    return out


# -------------------------------------------------------------- conditions

def _cmp(op, a, b):
    if op == "==":
        return a == b
    if op == "!=":
        return a != b
    if op == ">":
        return a > b
    if op == ">=":
        return a >= b
    if op == "<":
        return a < b
    if op == "<=":
        return a <= b
    if op == "in":
        return a in b
    if op == "not_in":
        return a not in b
    if op in ("contains", "not_contains"):
        try:
            hit = b in a
        except TypeError as e:
            raise ExprError(f"'{op}' needs a container/string element, "
                            f"got {a!r}: {e}") from None
        return hit if op == "contains" else not hit
    raise ExprError(f"unknown condition op {op!r}")


def _cond_ok(row, cond, env):
    if not isinstance(cond, (list, tuple)) or len(cond) not in (3, 4):
        raise ExprError(f"condition must be [field, op, value(, coerce)]: {cond!r}")
    field, op, value = cond[0], cond[1], cond[2]
    coerce = cond[3] if len(cond) == 4 else None
    if not isinstance(row, dict):
        raise ExprError(f"count-where row is {type(row).__name__}, not a dict")
    if field not in row:
        raise ExprError(f"count-where field {field!r} absent from row "
                        f"(have {sorted(row)[:8]}...)")
    v = row[field]
    if coerce == "int":
        v = int(v)
    elif coerce == "float":
        v = float(v)
    elif coerce is not None:
        raise ExprError(f"unknown coerce {coerce!r} (int/float)")
    if isinstance(value, dict) and set(value) == {"$ref"}:
        ref = value["$ref"]
        if ":" not in ref:
            raise ExprError(f"$ref must be 'alias:/pointer', got {ref!r}")
        alias, ptr = ref.split(":", 1)
        value = resolve(env.get(alias), ptr)
    return _cmp(op, v, value)


# -------------------------------------------------------------- functions

def _make_functions(env):
    def f(alias, pointer):
        return resolve(env.get(alias), pointer)

    def vals(alias, pattern):
        return resolve_pattern(env.get(alias), pattern)

    def rows(alias, pointer=""):
        r = resolve(env.get(alias), pointer)
        if not isinstance(r, list):
            raise ExprError(f"rows() needs a list at {alias}:{pointer}, "
                            f"got {type(r).__name__}")
        return r

    def cw(alias, pointer, conds):
        if isinstance(conds, str):
            try:
                conds = json.loads(conds)
            except json.JSONDecodeError as e:
                raise ExprError(f"cw conditions not valid JSON: {e}") from None
        if not isinstance(conds, list):
            raise ExprError("cw conditions must be a list of conditions")
        return sum(1 for row in rows(alias, pointer)
                   if all(_cond_ok(row, c, env) for c in conds))

    def where(xs, op, v):
        return [x for x in xs if _cmp(op, x, v)]

    def split(s, sep, i):
        try:
            return str(s).split(sep)[i]
        except IndexError:
            raise ExprError(f"split index {i} out of range for {s!r}") from None

    def jsondump(x):
        return json.dumps(x, sort_keys=True)

    def join(xs, sep):
        return sep.join(str(x) for x in xs)

    def joinrows(alias, pointer, template, sep):
        out = []
        for r in rows(alias, pointer):
            if not isinstance(r, dict):
                raise ExprError("joinrows rows must be dicts")
            try:
                out.append(template.format(**r))
            except KeyError as e:
                raise ExprError(f"joinrows template key {e} absent") from None
        return sep.join(out)

    def zipjoin(xs, ys, template, sep):
        if not isinstance(xs, list) or not isinstance(ys, list):
            raise ExprError(f"zipjoin needs two lists, got "
                            f"{type(xs).__name__}/{type(ys).__name__}")
        if len(xs) != len(ys):
            raise ExprError(f"zipjoin length mismatch: "
                            f"{len(xs)} vs {len(ys)} — parallel list fields "
                            "must pair row-for-row")
        try:
            return sep.join(template.format(a=x, b=y)
                            for x, y in zip(xs, ys))
        except (KeyError, IndexError, ValueError) as e:
            raise ExprError(f"zipjoin template failed: {e}") from None

    def frac(x):
        if isinstance(x, Fraction):
            return x
        if isinstance(x, bool) or not isinstance(x, (str, int, float)):
            raise ExprError(f"frac() takes a string/int/float, got {x!r}")
        try:
            return Fraction(str(x))
        except (ValueError, ZeroDivisionError) as e:
            raise ExprError(f"frac() cannot parse {x!r}: {e}") from None

    return {
        "f": f, "vals": vals, "rows": rows, "cw": cw, "where": where,
        "split": split, "jsondump": jsondump, "join": join,
        "joinrows": joinrows, "zipjoin": zipjoin, "frac": frac,
        "len": len, "sum": sum, "min": min, "max": max, "abs": abs,
        "int": int, "float": float, "str": str, "round": round, "bool": bool,
    }


# -------------------------------------------------------------- evaluator

_BINOPS = {
    ast.Add: lambda a, b: a + b,
    ast.Sub: lambda a, b: a - b,
    ast.Mult: lambda a, b: a * b,
    ast.Div: lambda a, b: a / b,
    ast.FloorDiv: lambda a, b: a // b,
    ast.Mod: lambda a, b: a % b,
    ast.Pow: lambda a, b: a ** b,
}
_CMPOPS = {
    ast.Eq: "==", ast.NotEq: "!=", ast.Lt: "<", ast.LtE: "<=",
    ast.Gt: ">", ast.GtE: ">=", ast.In: "in", ast.NotIn: "not_in",
}


def evaluate(expr, env, variables=None):
    """Evaluate one field expression against the receipt Env.

    `variables` maps names (a claim's `lets`) to already-computed values.
    """
    if not isinstance(expr, str):
        raise ExprError(f"expression must be a string, got {expr!r}")
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as e:
        raise ExprError(f"syntax error in expression {expr!r}: {e}") from None
    funcs = _make_functions(env)
    names = dict(variables or {})

    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant):
            if node.value is None or isinstance(node.value, (str, int, float, bool)):
                return node.value
            raise ExprError(f"disallowed constant {node.value!r}")
        if isinstance(node, ast.Name):
            if node.id in names:
                return names[node.id]
            raise ExprError(f"unknown name {node.id!r} (not in lets)")
        if isinstance(node, ast.Call):
            if node.keywords:
                raise ExprError("keyword arguments are not allowed")
            if not isinstance(node.func, ast.Name) or node.func.id not in funcs:
                raise ExprError("only whitelisted function calls are allowed")
            return funcs[node.func.id](*[ev(a) for a in node.args])
        if isinstance(node, ast.BinOp):
            op = _BINOPS.get(type(node.op))
            if op is None:
                raise ExprError(f"disallowed operator {type(node.op).__name__}")
            return op(ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return -ev(node.operand)
            if isinstance(node.op, ast.UAdd):
                return +ev(node.operand)
            if isinstance(node.op, ast.Not):
                return not ev(node.operand)
            raise ExprError(f"disallowed unary {type(node.op).__name__}")
        if isinstance(node, ast.Compare):
            left = ev(node.left)
            for op, comp in zip(node.ops, node.comparators):
                sym = _CMPOPS.get(type(op))
                if sym is None:
                    raise ExprError(f"disallowed comparison {type(op).__name__}")
                right = ev(comp)
                if not _cmp(sym, left, right):
                    return False
                left = right
            return True
        if isinstance(node, ast.BoolOp):
            if isinstance(node.op, ast.And):
                r = True
                for v in node.values:
                    r = ev(v)
                    if not r:
                        return r
                return r
            r = False
            for v in node.values:
                r = ev(v)
                if r:
                    return r
            return r
        if isinstance(node, (ast.List, ast.Tuple)):
            return [ev(e) for e in node.elts]
        if isinstance(node, ast.Dict):
            if any(k is None for k in node.keys):
                raise ExprError("dict unpacking (**) is not allowed")
            return {ev(k): ev(v) for k, v in zip(node.keys, node.values)}
        raise ExprError(f"disallowed syntax {type(node).__name__}")

    return ev(tree)
