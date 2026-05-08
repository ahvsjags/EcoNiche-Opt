from econiche_opt import __version__
from econiche_opt.cli import build_parser


def test_import_version_and_cli_commands():
    assert __version__
    parser = build_parser()
    help_text = parser.format_help()
    assert "validate-registry" in help_text
    assert "validate-project" in help_text
    assert "run-response-composite" in help_text
