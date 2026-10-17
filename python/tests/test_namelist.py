import pytest

from umwm import namelist

SAMPLE = """
! leading comment
&DOMAIN
  isGlobal     = .false.               ! Global (.T.) or regional (.F.)
  mm           = 21                    ! Domain size in x
  fmin         = 0.0313
  startTimeStr = '2026-01-01 00:00:00' ! not a ! comment
  dtg          = 3600
  flags        = .T. .F. T F
  restart      = .false.
/
&STOKES
  depths = 0.1 0.2, 0.3
           1.0e2 1d1
/
&PAIRS
  a = 1, b = 'two', c = "th""ree"
/
"""


def test_loads_groups_and_values():
    groups = namelist.loads(SAMPLE)
    assert list(groups) == ["domain", "stokes", "pairs"]
    domain = groups["domain"]
    assert domain["isglobal"] is False
    assert domain["mm"] == 21 and isinstance(domain["mm"], int)
    assert domain["fmin"] == pytest.approx(0.0313)
    assert domain["starttimestr"] == "2026-01-01 00:00:00"
    assert domain["dtg"] == 3600
    assert domain["flags"] == [True, False, True, False]
    assert groups["stokes"]["depths"] == [0.1, 0.2, 0.3, 100.0, 10.0]
    assert groups["pairs"] == {"a": 1, "b": "two", "c": 'th"ree'}


def test_loads_repository_namelist(repo_root):
    groups = namelist.load(repo_root / "namelists" / "main.nml")
    assert list(groups) == ["domain", "physics", "grid", "forcing", "forcing_constant", "output", "stokes"]
    assert groups["domain"]["om"] == 37
    assert groups["physics"]["nu_air"] == pytest.approx(1.56e-5)
    assert groups["output"]["stokes"] is True
    assert len(groups["stokes"]["depths"]) == 42


def test_group_terminated_on_value_line():
    groups = namelist.loads("&A\n x = 1 /\n")
    assert groups == {"a": {"x": 1}}


def test_unterminated_group_raises():
    with pytest.raises(namelist.NamelistError):
        namelist.loads("&A\n x = 1\n")


def test_value_without_name_raises():
    with pytest.raises(namelist.NamelistError):
        namelist.loads("&A\n 1\n/\n")


def test_dumps_round_trip():
    groups = namelist.loads(SAMPLE)
    text = namelist.dumps(groups)
    assert namelist.loads(text) == groups
    assert text.startswith("&DOMAIN\n")
    assert "  starttimestr   = '2026-01-01 00:00:00'" in text


@pytest.mark.parametrize("value, expected", [
    (True, ".true."), (False, ".false."), (3, "3"), (2.0, "2.0"), (1.56e-05, "1.56e-05"),
    ("it's", "'it''s'"), ([0.1, 1.0], "0.1 1.0"),
])
def test_format_value(value, expected):
    assert namelist.format_value(value) == expected


def test_dumps_with_comments():
    text = namelist.dumps({"g": {"x": 1}}, comments={("g", "x"): "a value"})
    assert "! a value" in text
    assert namelist.loads(text) == {"g": {"x": 1}}
