# Workspace Path Conventions

- 所有工作代码默认放在 `projectmd/`。
- 原始数据默认从 `projectfile/` 读取，视为只读目录。
- 论文和参考资料放在 `read/`，不要向该目录写入运行结果。
- 所有运行产物都写入 `output/`，包括图片、导出表格、h5ad、rds、中间结果和最终结果。
- 不要硬编码 `/mnt/...`、`C:/...` 等绝对路径。
- Python 代码优先使用 `workspace_paths.py` 中的 `raw_data_path()` 和 `output_path()`。
- R 代码优先 `source("workspace_paths.R")`，并使用 `raw_data_path()` 和 `output_path()`。
- 如果没有特殊说明，读取路径默认指向 `projectfile/`，写入路径默认指向 `output/`。