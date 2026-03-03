import sys
from pathlib import Path

from click.testing import CliRunner

project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root / "src"))

from lotus.cli import cli
from lotus.commands.inject import format_code_string, format_file, format_paths


def test_format_code_string_adds_logged_and_import():
    source = """def foo(adata):\n    return adata\n"""
    output = format_code_string(source)

    assert "from lotus.lineagetracker import logged" in output
    assert "@logged" in output
    assert "def foo(adata):" in output


def test_format_code_string_respects_prefix():
    source = """def keep_a(adata):\n    return adata\n\ndef skip_b(adata):\n    return adata\n"""
    output = format_code_string(source, target_prefix="keep_")

    assert "@logged\ndef keep_a" in output
    assert "@logged\ndef skip_b" not in output


def test_format_code_string_is_idempotent_and_no_duplicate_import():
    source = """from lotus.lineagetracker import logged\n\n@logged\ndef foo(adata):\n    return adata\n"""
    output = format_code_string(source)

    assert output.count("from lotus.lineagetracker import logged") == 1
    assert output.count("@logged") == 1


def test_format_file_and_format_paths(tmp_path):
    target_file = tmp_path / "mod.py"
    target_file.write_text("def foo(adata):\n    return adata\n", encoding="utf-8")
    changed = format_file(target_file)
    assert changed is True

    nested_dir = tmp_path / "pkg"
    nested_dir.mkdir()
    nested_file = nested_dir / "sub.py"
    nested_file.write_text("def bar(adata):\n    return adata\n", encoding="utf-8")

    ignored_dir = tmp_path / ".venv"
    ignored_dir.mkdir()
    ignored_file = ignored_dir / "ignored.py"
    ignored_file.write_text("def skip(adata):\n    return adata\n", encoding="utf-8")

    stats = format_paths(tmp_path)
    assert stats["scanned_files"] == 2
    assert stats["changed_files"] == 1
    assert "@logged" in nested_file.read_text(encoding="utf-8")
    assert "@logged" not in ignored_file.read_text(encoding="utf-8")


def test_cli_group_and_inject_command(tmp_path):
    file_path = tmp_path / "target.py"
    file_path.write_text("def keep_me(adata):\n    return adata\n", encoding="utf-8")

    runner = CliRunner()
    help_result = runner.invoke(cli, ["--help"])
    assert help_result.exit_code == 0
    assert "inject" in help_result.output

    result = runner.invoke(cli, ["inject", str(file_path), "--prefix", "keep_"])
    assert result.exit_code == 0
    updated = file_path.read_text(encoding="utf-8")
    assert "@logged" in updated
    assert "from lotus.lineagetracker import logged" in updated
