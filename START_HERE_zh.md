# 仓库初始化说明

本包是二维热传导项目的基础骨架，包含说明、构建配置、实验参数和文档框架。当前没有实现求解器，也没有虚构任何测试或性能结果。

## GitHub 的基本信息

- 建议仓库名：`parallel-heat2d`。
- 仓库 description：复制根目录 `REPOSITORY_DESCRIPTION.txt` 中的一行英文。
- README：使用本包中的 `README.md`。
- 可选 topics：`cpp`、`mpi`、`openmp`、`parallel-computing`、`heat-equation`、`hpc`。

description 和 topics 要在 GitHub 仓库页面中手动设置；提交同名文本文件不会自动修改 GitHub 的这些字段。

## 如何放进新仓库

1. 在本地克隆你刚建立的 GitHub 仓库。
2. 解压本包。
3. 把解压后 `parallel-heat2d` 文件夹里的内容复制到本地仓库根目录。根目录应直接看到 `README.md` 和 `CMakeLists.txt`，不要再多套一层文件夹。
4. 同时复制 `.gitignore`、`.gitattributes`、`.editorconfig`。
5. 如果 GitHub 初始化时已经生成 README 或 `.gitignore`，先比较内容，再合并或替换。
6. 在 VS Code 的 Source Control 中检查文件，commit 后 push。

也可以在本地仓库根目录执行下面的命令。先阅读 `git status`，确认即将提交的内容正是这次的基础文件。

```bash
git status
git add README.md START_HERE_zh.md REPOSITORY_DESCRIPTION.txt CONTRIBUTING.md CMakeLists.txt .gitignore .gitattributes .editorconfig requirements.txt config experiments src scripts hpc results docs report
git commit -m "Initialize heat diffusion project scaffold"
git push
```

这些命令不会创建远程仓库；需要先完成克隆或远程关联。如果尚未设置上游分支，Git 会提示对应的首次 push 命令。

## 当前哪些东西能用

| 内容 | 状态 |
| --- | --- |
| 英文 README、仓库 description | 已提供 |
| 默认参数和实验配置 | 已提供，属于后续实现的接口约定 |
| CMake 骨架配置 | 已提供，默认不编译求解器 |
| 串行及 MPI/OpenMP 求解器 | 尚未实现 |
| 批量实验、结果检查及绘图工具 | 尚未实现，目录里有工具职责说明 |
| HPC 作业模板 | 已提供；需要按实际环境设置模块、资源和启动命令 |
| 实验结果及报告结论 | 尚未生成 |

有 CMake 和 C++ 编译器后，可以检查基础工程配置：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

当前不会生成 `heat_serial` 或 `heat_parallel`。README 里的求解器运行命令是后续要实现的目标接口。

## 文件的主要用途

- `config/defaults.json`：默认物理条件、网格尺寸、步数和输出设置。
- `experiments/*.json`：正确性、强扩展、弱扩展和策略比较的计划配置。
- `src/common/README.md`：公共数据结构和更新公式怎么组织。
- `src/serial/README.md`、`src/parallel/README.md`：两种求解器的实现边界。
- `docs/design.md`：需要统一的符号、索引、边界和通信规则。
- `docs/experiment_plan.md`：如何验证、计时和计算加速比。
- `docs/ai_usage.md`：已经发生的 AI 辅助及后续核查记录。
- `report/report.md`：最终报告框架，待补充真实结果。

## 下一步的实现顺序

先完成公共网格和参数，再实现串行求解器及正确性检查。随后依次加入按行划分、MPI 边界通信、OpenMP、二维分块和通信计算重叠。

建议使用分支开发，并让同一批公共公式服务于串行和并行版本。共享公式不能代替独立的物理验证，因此配置中也保留了恒温场和正弦衰减用例。
