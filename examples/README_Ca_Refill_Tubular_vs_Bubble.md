# Ca²⁺ Refill Examples: Tubular vs Bubble ER Networks

This folder contains example input files and workflows for the **Ca²⁺ refill project**, comparing
a **tubular ER network** and a **bubble-based ER network** under the same cell-scale geometry.

---

## Network files (in `examples/`)

### Tubular ER network
- **`diamondR10_l1_sphere.net`**
  - Cell radius: **10 µm**
  - Perinuclear sheet radius: **5 µm**
  - ER geometry: **diamond tubular network**
  - Tube radius: **0.05 µm**

### Bubble ER network
- **`HCPdeg6_bubblesRpt6_ellpt3.net`**
  - Same cell radius (**10 µm**) and perinuclear sheet radius (**5 µm**)
  - ER geometry: **connected bubble network**
  - Bubble radius: **0.6 µm**
  - Connecting tube radius: **0.018 µm**

---

## Parameter files

- **Tubular network**
  - `param.diamondR10_nf_nuccon91.0.0`
- **Bubble network**
  - `param.HCP_nf_bub_10sheet_buffscl_R10.0.0`

These files specify diffusion constants, buffering parameters, reservoir geometry,
and which network nodes are treated as fixed (PM-contact) nodes.

---

## How to run

### 1. Compile the code
1. Make sure you have a **Fortran compiler** (e.g. `gfortran`) installed.
2. Go into the `source/` directory and run:
   ```bash
   make
   ```
   This will generate the executable:
   ```
   netmeshdynamicsFVM.exe
   ```

---

### 2. Run Ca²⁺ refill simulations
From the `examples/` directory:

```bash
../netmeshdynamicsFVM.exe diamondR10_nf_nuccon91.0.0 > stdout.diamondR10_nf_nuccon91.0.0
../netmeshdynamicsFVM.exe HCP_nf_bub_10sheet_buffscl_R10.0.0 > stdout.HCP_nf_bub_10sheet_buffscl_R10.0.0
```

Each run will produce several output files, including:
- `*.snap.txt`: contains the concentration of free calcium and total buffer protein at all mesh cells, at each snapshot time
- `*.mesh.txt`
- `*.out`: contains the total flux of calcium out of the network at each time point

---

## Post-processing and plotting

A MATLAB script is provided to parse and visualize the Ca²⁺ refill results:

- **`examples/plot_Ca_Refill.m`**

This script:
- Reads `*.snap.txt` and `*.mesh.txt` files
- Extracts ER Ca²⁺ concentration over time
- Produces plots comparing **tubular vs bubble networks**

---

## Summary

These examples demonstrate how **ER geometry**
(tubular vs bubble-based networks) affects Ca²⁺ refill dynamics under
otherwise identical cellular and biochemical conditions.
