try:
    from importlib.metadata import version as _get_version

    pyslim_version = _get_version("pyslim")
except Exception:
    pyslim_version = "unknown"

slim_file_version = "0.9"
# other file versions that require no modification
compatible_slim_file_versions = ["0.9"]
