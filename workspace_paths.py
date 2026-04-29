from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parent
READ_DIR = WORKSPACE_ROOT / "read"
RAW_DATA_DIR = WORKSPACE_ROOT / "projectfile"
OUTPUT_DIR = WORKSPACE_ROOT / "output"
CODE_DIR = WORKSPACE_ROOT / "projectmd"


def workspace_dirs() -> dict[str, Path]:
    return {
        "root": WORKSPACE_ROOT,
        "read": READ_DIR,
        "projectfile": RAW_DATA_DIR,
        "output": OUTPUT_DIR,
        "projectmd": CODE_DIR,
    }


def read_path(*parts: str) -> Path:
    return READ_DIR.joinpath(*parts)


def raw_data_path(*parts: str) -> Path:
    return RAW_DATA_DIR.joinpath(*parts)


def output_path(*parts: str) -> Path:
    path = OUTPUT_DIR.joinpath(*parts)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def code_path(*parts: str) -> Path:
    return CODE_DIR.joinpath(*parts)