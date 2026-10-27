# VS Code + WSL Ubuntu：安装、编译和运行

这次已经提供真实的 C++ 和 Python 实现。串行、MPI 按行划分、MPI 二维分块、OpenMP、阻塞通信及通信计算重叠均已包含。

本说明按 Windows + WSL Ubuntu 编写。除了第一步的 PowerShell 命令，后面的命令都在 **Ubuntu / VS Code 的 WSL 终端**运行。

## 1. 准备 WSL 和 VS Code

如果已经有 WSL Ubuntu，直接进入下一步。否则，在管理员 PowerShell 中执行：

```powershell
wsl --install -d Ubuntu
```

按提示重启，第一次打开 Ubuntu 时设置 Linux 用户名和密码。在 Windows 安装 VS Code，并安装微软的 WSL 扩展。打开项目后，按提示在 WSL 环境安装 C/C++ 和 CMake Tools。

在 Ubuntu 里执行 `code .` 会用 VS Code 打开当前目录。左下角应显示 `WSL: Ubuntu` 或对应发行版名称。

## 2. 放置仓库文件

建议把项目放在 Ubuntu 文件系统，例如 `~/projects/parallel-heat2d`。如果需要从 Windows 复制解压后的文件夹，可以先在 Ubuntu 运行：

```bash
mkdir -p ~/projects
cd ~/projects
explorer.exe .
```

资源管理器会打开这个 Linux 目录。把解压后的 `parallel-heat2d` 文件夹放进去；如果已经有本地 Git 仓库，将本包内容合并到那个仓库，保留原来的 `.git`。

进入项目根目录并打开 VS Code：

```bash
cd ~/projects/parallel-heat2d
code .
```

根目录应直接包含 `CMakeLists.txt`、`README.md`、`src` 和 `scripts`。复制时包含 `.vscode`、`.gitignore` 等文件。

## 3. 安装 Ubuntu 编译与运行环境

在 Ubuntu 或 VS Code 的 WSL 终端执行：

```bash
sudo apt update
sudo apt install -y build-essential cmake gdb openmpi-bin libopenmpi-dev python3-venv python3-pip
```

检查命令是否可用：

```bash
g++ --version
cmake --version
mpicxx --version
mpiexec --version
```

GCC 编译 C++，CMake 管理构建，Open MPI 提供 MPI 通信和多进程启动。GCC 的 OpenMP 支持会由 CMake 自动启用。

## 4. 准备 Python 工具环境

在仓库根目录执行：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

这里的 Python 只用于检查、实验调度和绘图，热传导求解在 C++ 程序中完成。以后打开新的终端，使用 Python 工具前重新执行 `source .venv/bin/activate`。

## 5. 编译项目

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DHEAT2D_BUILD_SOLVERS=ON
cmake --build build --parallel 2
```

这会生成：

- `build/heat_serial`：纯串行程序。
- `build/heat_parallel`：MPI + OpenMP 程序，通过参数切换并行策略。

如果以前配置过旧骨架，显式的 `HEAT2D_BUILD_SOLVERS=ON` 会覆盖缓存中的 OFF。修改 `.cpp` 或 `.hpp` 后，再运行第二条构建命令。

如果要练习断点调试，可以按 F5，选择 `Debug serial solver (WSL)`；配置会在 `build-debug` 中另建 Debug 版本。Ctrl+Shift+B 对应 Release 构建任务。当前 F5 配置用于串行调试，多进程验证用终端命令。

## 6. 运行串行版

```bash
./build/heat_serial --nx 128 --ny 128 --steps 2000 --case heated_top --output outputs/serial.csv
```

含义：128×128 个内部点，推进 2000 步，上边界加热，最后导出温度场。CSV 包含物理边界，所以是 130×130。终端会打印一行 JSON，包含实际时间步长、计算耗时等信息。

`nx`、`ny` 是空间网格大小，`steps` 是时间步数。`cfl` 默认 0.8，程序会据此选择稳定时间步长；实际物理时长等于 `steps × dt`。

## 7. 运行并行版

先试 2 个进程，每个进程 2 个线程：

```bash
OMP_NUM_THREADS=2 mpiexec --bind-to none -n 2 ./build/heat_parallel --nx 128 --ny 128 --steps 2000 --case heated_top --decomp 2d --exchange overlap --output outputs/parallel.csv
```

这条本地演示命令启动 4 个 CPU 工作线程。如果资源较少，可以把 `OMP_NUM_THREADS` 改成 1。这里 `--bind-to none` 是 Open MPI 的本地演示选项，不要直接拿这个配置宣称最佳性能。

可切换的参数：

| 参数 | 取值 | 含义 |
| --- | --- | --- |
| `--decomp` | `rows` | 按行切分 |
| `--decomp` | `2d` | 二维分块 |
| `--exchange` | `blocking` | 完成通信后再计算 |
| `--exchange` | `overlap` | 先发起通信，计算内部，再等待并计算边缘 |
| `--py`、`--px` | 正整数 | 可显式指定进程网格，乘积等于进程数 |

`OMP_NUM_THREADS=1` 可用于 MPI 单线程基准。需要把进程与线程都分配到真实可用 CPU 上，最终性能实验的绑定方式要根据机器配置确认。

## 8. 自动检查正确性

```bash
python scripts/check_correctness.py --quick
```

快速检查最多使用 2 个进程、每个进程 1 个线程，会检查手算的一步更新、恒温场、固定边界、串并行一致性、无效参数，以及正弦温度衰减的网格收敛。

有足够资源时再运行完整配置：

```bash
python scripts/check_correctness.py
python tests/test_integration.py
```

完整检查最多使用 6 个进程和每进程 2 个线程。扩展检查还覆盖很窄的网格、大消息、二维列通信和不合法的进程划分。输出记录在 `outputs/validation/summary.json` 和 `outputs/extended_validation/summary.json`。

## 9. 查看温度图

```bash
python scripts/plot_temperature.py outputs/serial.csv --output outputs/temperature.png --vmin 20 --vmax 100
```

在 VS Code 左侧点击 `outputs/temperature.png` 即可查看。`examples` 中也附了一张实际生成的示例图。

## 10. 检查批量实验流程

先用很小的配置跑通整个链路：

```bash
python scripts/run_benchmarks.py --config experiments/smoke_scaling.json
python scripts/plot_results.py outputs/smoke/metrics.csv --out-dir outputs/smoke/figures
```

它会产生原始计时 CSV、每次运行日志、完整配置，以及运行时间、加速比和效率图。这只是本地流程检查，不能代替 Atlas/Vanda 的最终实验。

完整实验分别使用：

- `experiments/strong_scaling.json`：固定问题，增加核数。
- `experiments/weak_scaling.json`：固定每个进程的局部网格及线程数，扩大问题。
- `experiments/strategies.json`：固定总核数，比较进程/线程组合、分区和通信策略。

运行大实验前，可先解析配置：

```bash
python scripts/run_benchmarks.py --config experiments/strong_scaling.json --dry-run
```

正式运行时，必须匹配可用 CPU、内存和进程绑定。README 提供了 Open MPI 的按每进程线程数绑定示例；HPC 队列、模块和启动方式仍需按实际平台填写。

## 11. 常见情况

- 找不到 `cmake` 或 `mpicxx`：确认第 3 步在 Ubuntu 内执行，并检查终端属于 WSL。
- 找不到 `numpy`：激活 `.venv`，再安装 `requirements.txt`。
- 只配置成功却没有程序：重新配置并显式传入 `-DHEAT2D_BUILD_SOLVERS=ON`，然后执行构建。
- MPI 提示 slots 不足：减少进程数或在正确资源分配内运行。本地正确性演示可以考虑 Open MPI 的 `--oversubscribe`，但不要用超额订阅的数据评价并行效率。
- 希望增加运行步数：修改 `--steps`，不要随意增大 `--dt`；程序会拒绝不稳定步长。
- 使用之前的 Git 仓库：合并源文件后检查 `git status`，再提交；本包不含 `.git`，不会提供你的远程地址或账户信息。

目前已完成本地 Linux 数值与脚本测试；这里没有实际操作 Windows 桌面的 WSL 安装，也没有提交 Atlas/Vanda 作业。实际验证范围见 `docs/verification.md`。
