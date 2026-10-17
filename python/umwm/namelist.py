"""Reading and writing of Fortran namelist files.

This module handles the subset of the namelist syntax used by the UMWM
configuration files: named groups ``&GROUP ... /`` holding ``name = value``
assignments, where a value is a logical, integer, real, or quoted string, or a
list of these separated by spaces or commas; ``!`` starts a comment. Names
are case-insensitive and are returned in lower case.
"""

from __future__ import annotations

import re
from collections import OrderedDict
from typing import Any

__all__ = ["loads", "load", "dumps", "dump", "format_value", "parse_scalar", "parse_values"]

_GROUP_START = re.compile(r"^\s*&\s*(\w+)\s*$")
_TOKEN = re.compile(
    r"""
    '(?:[^']|'')*'          # single-quoted string
    | "(?:[^"]|"")*"        # double-quoted string
    | [^\s,=]+              # bare value or name
    | =
    | ,
    """,
    re.VERBOSE,
)
_LOGICAL_TRUE = {".true.", ".t.", "t", "true"}
_LOGICAL_FALSE = {".false.", ".f.", "f", "false"}
_INTEGER = re.compile(r"^[+-]?\d+$")
_REAL = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")


class NamelistError(ValueError):
    """Raised for a namelist file that cannot be parsed."""


def _strip_comment(line: str) -> str:
    """Remove a trailing ``!`` comment, honoring quoted strings."""
    quote = None
    for i, char in enumerate(line):
        if quote:
            if char == quote:
                quote = None
        elif char in "'\"":
            quote = char
        elif char == "!":
            return line[:i]
    return line


def parse_scalar(token: str) -> Any:
    """Convert a single namelist value token to a Python object."""
    if len(token) >= 2 and token[0] == token[-1] and token[0] in "'\"":
        quote = token[0]
        return token[1:-1].replace(quote * 2, quote)
    lowered = token.lower()
    if lowered in _LOGICAL_TRUE:
        return True
    if lowered in _LOGICAL_FALSE:
        return False
    if _INTEGER.match(token):
        return int(token)
    if _REAL.match(token):
        return float(token.lower().replace("d", "e"))
    return token


def parse_values(text: str) -> list[Any]:
    """Convert a namelist value list (values separated by spaces or commas,
    strings quoted) to a list of Python objects."""
    return [parse_scalar(token) for token in _TOKEN.findall(text) if token not in (",", "=")]


def loads(text: str) -> "OrderedDict[str, OrderedDict[str, Any]]":
    """Parse namelist text into an ordered mapping of groups to variables.

    A variable with several values is returned as a list; a variable with a
    single value is returned as a scalar.
    """
    groups: OrderedDict[str, OrderedDict[str, Any]] = OrderedDict()
    group: OrderedDict[str, Any] | None = None
    group_name = None
    current = None
    values: list[Any] = []

    def finish_variable():
        nonlocal current, values
        if current is not None:
            group[current] = values[0] if len(values) == 1 else list(values)
        current = None
        values = []

    for number, raw in enumerate(text.splitlines(), start=1):
        line = _strip_comment(raw).strip()
        if not line:
            continue
        if group is None:
            match = _GROUP_START.match(line)
            if not match:
                if line.startswith("&"):
                    raise NamelistError(f"line {number}: malformed group header: {raw!r}")
                continue  # text outside groups is ignored, as Fortran does
            group_name = match.group(1).lower()
            group = OrderedDict()
            continue
        if line == "/" or line == "&end" or line.lower() == "&end":
            finish_variable()
            groups[group_name] = group
            group = None
            continue
        tokens = _TOKEN.findall(line)
        if line.endswith("/"):
            tokens = _TOKEN.findall(line[:-1])
            end_of_group = True
        else:
            end_of_group = False
        index = 0
        while index < len(tokens):
            token = tokens[index]
            if token == ",":
                index += 1
                continue
            if index + 1 < len(tokens) and tokens[index + 1] == "=":
                finish_variable()
                current = token.lower()
                index += 2
                continue
            if current is None:
                raise NamelistError(f"line {number}: value without a variable name: {raw!r}")
            values.append(parse_scalar(token))
            index += 1
        if end_of_group:
            finish_variable()
            groups[group_name] = group
            group = None

    if group is not None:
        raise NamelistError(f"group &{group_name} is not terminated by '/'")
    return groups


def load(path) -> "OrderedDict[str, OrderedDict[str, Any]]":
    """Parse the namelist file at ``path``."""
    with open(path, encoding="utf-8") as stream:
        return loads(stream.read())


def format_value(value: Any) -> str:
    """Format a Python object as a namelist value."""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, (list, tuple)):
        if not value:
            raise NamelistError("an empty list cannot be written to a namelist")
        return " ".join(format_value(item) for item in value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        text = repr(float(value))
        return text if ("e" in text or "." in text or "n" in text) else text + "."
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    raise NamelistError(f"cannot write a value of type {type(value).__name__} to a namelist")


def dumps(groups, comments=None, width=14) -> str:
    """Format a mapping of groups to variables as namelist text.

    ``comments`` optionally maps ``(group, variable)`` to a trailing comment.
    """
    comments = comments or {}
    lines = []
    for group_name, variables in groups.items():
        lines.append(f"&{group_name.upper()}")
        for name, value in variables.items():
            entry = f"  {name:<{width}} = {format_value(value)}"
            comment = comments.get((group_name, name))
            if comment:
                entry = f"{entry:<44} ! {comment}"
            lines.append(entry)
        lines.append("/")
    return "\n".join(lines) + "\n"


def dump(groups, path, comments=None) -> None:
    """Write a mapping of groups to variables to the namelist file at ``path``."""
    with open(path, "w", encoding="utf-8") as stream:
        stream.write(dumps(groups, comments=comments))
