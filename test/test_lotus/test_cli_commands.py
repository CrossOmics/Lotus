import sys
from pathlib import Path

from click.testing import CliRunner

project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root / "src") not in sys.path:
    sys.path.insert(0, str(project_root / "src"))

from lotus.cli import cli


def _run_inject(tmp_path, source_code: str, prefix: str = "") -> str:
    target = tmp_path / "target.py"
    target.write_text(source_code, encoding="utf-8")
    runner = CliRunner()
    args = ["inject", str(target)]
    if prefix:
        args.extend(["--prefix", prefix])
    result = runner.invoke(cli, args)
    assert result.exit_code == 0, result.output
    return target.read_text(encoding="utf-8")


def test_cli_inject_top_level_functions(tmp_path):
    source = (
        "def preprocess(adata):\n"
        "    return adata\n\n"
        "def summarize(adata):\n"
        "    return adata\n"
    )
    updated = _run_inject(tmp_path, source)

    assert "from lotus.lineagetracker import logged" in updated
    assert "@logged\ndef preprocess" in updated
    assert "@logged\ndef summarize" in updated


def test_cli_inject_nested_functions(tmp_path):
    source = (
        "def outer(adata):\n"
        "    def inner(x):\n"
        "        return x\n"
        "    return inner(adata)\n"
    )
    updated = _run_inject(tmp_path, source)

    assert "@logged\ndef outer" in updated
    assert "    @logged\n    def inner" not in updated
    assert "return inner(adata)" in updated


def test_cli_inject_class_fields_and_methods(tmp_path):
    source = (
        "class Processor:\n"
        "    version = 1\n\n"
        "    def __init__(self, adata):\n"
        "        self.adata = adata\n\n"
        "    def run(self):\n"
        "        return self.adata\n"
    )
    updated = _run_inject(tmp_path, source)

    assert "class Processor:" in updated
    assert "version = 1" in updated
    assert "def __init__" in updated
    assert "    @logged\n    def run" in updated


def test_cli_skips_virtual_env_directory(tmp_path):
    # Mark directory as a virtual env by creating pyvenv.cfg
    venv_dir = tmp_path / "myenv"
    venv_dir.mkdir()
    (venv_dir / "pyvenv.cfg").write_text("home = /usr/bin/python", encoding="utf-8")

    target = venv_dir / "mod.py"
    original = "def foo(adata):\n    return adata\n"
    target.write_text(original, encoding="utf-8")

    runner = CliRunner()
    result = runner.invoke(cli, ["inject", str(venv_dir)])

    assert result.exit_code == 0, result.output
    assert "Scanned 0 file(s), changed 0 file(s)." in result.output
    assert target.read_text(encoding="utf-8") == original
