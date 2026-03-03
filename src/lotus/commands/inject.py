"""Inject subcommand for adding @logged decorators."""

from __future__ import annotations

from pathlib import Path
import click

from typing import Iterable

import libcst as cst
import libcst.matchers as m

_DEFAULT_EXCLUDED_NAME = {
    ".git",
    ".venv",
    "__pycache__",
    "build",
    "dist",
    "__init__.py"
}


# ---- Helper Functions ----
def _is_logged_decorator(decorator: cst.Decorator) -> bool:
    return m.matches(
        decorator.decorator,
        m.OneOf(
            m.Name("logged"),
            m.Attribute(attr=m.Name("logged")),
        ),
    )

def is_virtual_env(path: Path) -> bool:
    # A virtual environment must be a directory
    if not path.is_dir():
        return False
        
    # Fingerprint 1: pyvenv.cfg (Standard for Python 3.3+ venv and modern virtualenv)
    if (path / "pyvenv.cfg").is_file():
        return True
        
    # Fingerprint 2: bin/activate (macOS/Linux) or Scripts/activate (Windows)
    if (path / "bin" / "activate").is_file() or (path / "Scripts" / "activate").is_file():
        return True
        
    return False

def _statement_has_logged_import(stmt: cst.CSTNode) -> bool:
    if not m.matches(stmt, m.SimpleStatementLine(body=[m.ImportFrom()])):
        return False
    line = stmt
    import_from = line.body[0]
    if not isinstance(import_from, cst.ImportFrom):
        return False
    names = import_from.names
    if isinstance(names, cst.ImportStar):
        return True
    for alias in names:
        if isinstance(alias, cst.ImportAlias) and m.matches(alias.name, m.Name("logged")):
            return True
    return False


def _module_has_logged_import(module: cst.Module) -> bool:
    return any(_statement_has_logged_import(stmt) for stmt in module.body)


def _is_import_stmt(stmt: cst.CSTNode) -> bool:
    return m.matches(
        stmt,
        m.SimpleStatementLine(
            body=[
                m.OneOf(
                    m.Import(),
                    m.ImportFrom(),
                )
            ]
        ),
    )


def _is_future_import(stmt: cst.CSTNode) -> bool:
    return m.matches(
        stmt,
        m.SimpleStatementLine(
            body=[
                m.ImportFrom(
                    module=m.Name("__future__"),
                )
            ]
        ),
    )


def _logged_import_statement() -> cst.SimpleStatementLine:
    return cst.SimpleStatementLine(
        body=[
            cst.ImportFrom(
                module=cst.Attribute(
                    value=cst.Name("lotus"),
                    attr=cst.Name("lineagetracker"),
                ),
                names=[cst.ImportAlias(name=cst.Name("logged"))],
            )
        ]
    )


def _insert_logged_import(module: cst.Module) -> cst.Module:
    # Insert the import after module docstring / __future__ imports / import block.
    body = list(module.body)
    insert_at = 0

    if body and m.matches(
        body[0],
        m.SimpleStatementLine(body=[m.Expr(value=m.SimpleString())]),
    ):
        insert_at = 1

    while insert_at < len(body) and _is_future_import(body[insert_at]):
        insert_at += 1

    while insert_at < len(body) and _is_import_stmt(body[insert_at]):
        insert_at += 1

    body.insert(insert_at, _logged_import_statement())
    return module.with_changes(body=body)


class LoggedDecoratorAdder(cst.CSTTransformer):
    """Add ``@logged`` to matching functions that don't already have it."""

    def __init__(self, target_prefix: str = ""):
        self.target_prefix = target_prefix
        self.changed = False

    def leave_FunctionDef(
        self,
        original_node: cst.FunctionDef,
        updated_node: cst.FunctionDef,
    ) -> cst.FunctionDef:
        # Exclude magic methods (e.g., __call__()) and methods don't start with target prefix (if there is)
        if (original_node.name.value.startswith("__") and original_node.name.value.endswith("__")) or \
            (self.target_prefix and not original_node.name.value.startswith(self.target_prefix)):
            return updated_node

        # If already contains @logged, skip the current function
        if any(_is_logged_decorator(decorator) for decorator in original_node.decorators):
            return updated_node

        self.changed = True
        new_decorator = cst.Decorator(decorator=cst.Name("logged"))
        return updated_node.with_changes(decorators=(new_decorator, *updated_node.decorators))


def format_code_string(source_code: str, target_prefix: str = "") -> str:
    """Return code with @logged injected, preserving formatting via LibCST."""
    module = cst.parse_module(source_code)
    transformer = LoggedDecoratorAdder(target_prefix=target_prefix)
    modified_module = module.visit(transformer)

    if transformer.changed and not _module_has_logged_import(modified_module):
        modified_module = _insert_logged_import(modified_module)

    return modified_module.code


def format_file(file_path: str | Path, target_prefix: str = "") -> bool:
    """Format one file in-place. Returns True only when file content changed."""
    path = Path(file_path)
    source_code = path.read_text(encoding="utf-8")
    modified_code = format_code_string(source_code, target_prefix)

    if source_code == modified_code:
        return False

    path.write_text(modified_code, encoding="utf-8")
    return True


def _iter_python_files(path: Path) -> Iterable[Path]:
    venv_cache: dict[Path, bool] = {}

    def _is_venv_cached(p: Path) -> bool:
        resolved = p.resolve()
        if resolved not in venv_cache:
            venv_cache[resolved] = is_virtual_env(resolved)
        return venv_cache[resolved]

    if path.is_file():
        if (
            path.suffix == ".py"
            and path.name not in _DEFAULT_EXCLUDED_NAME
            and not _is_venv_cached(path.parent)
        ):
            yield path
        return
    
    # Skip immediately when the target directory itself is a virtual env.
    if _is_venv_cached(path):
        return

    for candidate in path.rglob("*.py"):
        if any(part in _DEFAULT_EXCLUDED_NAME for part in candidate.parts):
            continue
        # Also skip any file nested under a detected virtual env directory.
        if any(_is_venv_cached(parent) for parent in candidate.parents):
            continue
        yield candidate


def format_paths(path: str | Path = ".", target_prefix: str = "") -> dict[str, int]:
    """Format all target Python files and return scan/change statistics."""
    target = Path(path)
    files = list(_iter_python_files(target))
    changed = 0
    for file_path in files:
        if format_file(file_path, target_prefix=target_prefix):
            changed += 1
    return {
        "scanned_files": len(files),
        "changed_files": changed,
    }

@click.command()
@click.argument(
    "path",
    required=False,
    default=".",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--prefix",
    default="",
    help="Only add @logged to functions whose names start with this prefix.",
)
def inject_cmd(path: Path, prefix: str) -> None:
    """Inject @logged decorators into a file or directory tree."""
    stats = format_paths(path=path, target_prefix=prefix)
    click.echo(
        f"Scanned {stats['scanned_files']} file(s), changed {stats['changed_files']} file(s)."
    )
