from mpro_docking.constants import EXHAUSTIVENESS, SEED
from mpro_docking.vina_io import format_vina_config, parse_score_only, parse_vina_modes, torsdof

LOG = """
Using random seed: 42
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1         -7.321          0          0
   2         -6.100      2.100      3.400
"""

PDBQT = """
REMARK VINA RESULT:      -7.321      0.000      0.000
REMARK VINA RESULT:      -6.100      2.100      3.400
TORSDOF 6
"""


def test_log_and_pdbqt_scores():
    log_modes = parse_vina_modes(LOG)
    pose_modes = parse_vina_modes(PDBQT)
    assert [mode.score for mode in log_modes] == [-7.321, -6.1]
    assert [mode.score for mode in pose_modes] == [-7.321, -6.1]
    assert pose_modes[0].rmsd_lb == 0.0
    assert torsdof(PDBQT) == 6


def test_score_only_line():
    text = "Estimated Free Energy of Binding    : -4.250 (kcal/mol)\n"
    assert parse_score_only(text) == -4.25


def test_config_records_the_protocol():
    text = format_vina_config(
        "structures/receptor/6LU7_receptor.pdbqt",
        (-10.729, 12.418, 68.816),
        (30.0, 30.0, 30.0),
    )
    assert "receptor = structures/receptor/6LU7_receptor.pdbqt" in text
    assert "center_x = -10.729" in text
    assert "center_y = 12.418" in text
    assert "center_z = 68.816" in text
    assert "size_x = 30" in text
    assert f"exhaustiveness = {EXHAUSTIVENESS}" in text
    assert f"seed = {SEED}" in text
