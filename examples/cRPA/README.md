# NiO cRPA 接口算例

本算例读取 `crpa/space.json`、`crpa/projectors.dat` 和 `crpa/UJ.dat`，
运行 NiO 的单杂质 DMFT。projector 由 NiO 原 NAO projector 导出；
相互作用张量由脚本生成，不是 cRPA 屏蔽计算的结果。

- 5 个 Ni-d 轨道，666 个 k 点，DFT band 编号 8–12（从 0 编号）。
- 温度 1160 K，投影窗口 −2 到 1 eV。
- U=4 eV、J=1 eV，density_density 相互作用，FLL 双计数。
- 固定 DFT，30 轮 DMFT，Rutgers CT-HYB 每条链 1,000,000 步。

`crpa/` 包含空间定义、projector、空间轨道张量和转换后的 U/J。
`dft/outputs_to_DMFT/` 包含能带、k 点权重、原子和对称性信息。
原子信息文件中的 U=1、J=0.1 是占位值，运行时由 `crpa/UJ.dat` 覆盖。
此路径不需要 AO 波函数或重叠矩阵文件。
projector 格式见[输入输出文档](../../../docs/input-output-files.md)。

编译 DMFT 并准备 MPI、Rutgers 求解器所需的运行环境后，在本目录执行：

```bash
python generate_Uijkl.py
python ../../../scripts/Uijkl_to_UJ.py
mpirun -np 4 ../../../bin/DFTDMFT > dmft.log 2>&1
python check_result.py
```

Python 需要 NumPy。输入文件已包含生成的张量和 U/J，因此前两步仅用于重新生成它们。
`generate_Uijkl.py` 使用空间轨道约定 `<ij|W|kl>`，写出完整五轨道张量：
同轨道 U=4，轨道间直接项 U−2J=2，交换与对跃迁项 J=1（单位 eV）。
