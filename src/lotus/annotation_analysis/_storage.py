from pathlib import Path


ANNOTATION_SCHEMA_VERSION = "1.0"


def sanitize_storage_key(value: str) -> str:
    return Path(value).stem.replace(" ", "_").replace("-", "_")


def build_obs_column_name(namespace: str, key: str, field: str) -> str:
    return f"{namespace}_{sanitize_storage_key(key)}_{field}"
