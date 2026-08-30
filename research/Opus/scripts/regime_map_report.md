# Regime Map Report: Mie-Gruneisen vs Ideal Gas EOS

This report maps the solution regimes of the 1-D elastoplastic piston problem as a function of the yield-to-shear ratio $Y_0/G$ and the piston velocity $v_{\text{piston}}$, for both the **original** (single total-energy) and **energy-split** modes, using two equations of state:

1. **Mie-Gruneisen** (MG): $P = P_H(\mu) + \Gamma_0 \rho (e - e_H(\mu))$ with Hugoniot reference curve from the linear $U_s = C_0 + s\,u_p$ fit.
2. **Ideal gas** (IG): $P = (\gamma - 1)\,\rho\,e$ with $\gamma = \Gamma_0 + 1 = 3$.

## Material parameters (CGS)

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Reference density | $\rho_0$ | 2.79 g/cm$^3$ |
| Bulk sound speed | $C_0$ | 5.330e+05 cm/s |
| Hugoniot slope | $s$ | 1.34 |
| Gruneisen parameter | $\Gamma_0$ | 2.0 |
| Shear modulus | $G$ | 2.860e+11 dyn/cm$^2$ |
| Yield stress | $Y_0$ | variable ($Y_0/G = 0.005$--$0.80$) |

## Solution regimes

The plastic-shock speed $U_s$ is determined by finding roots of

$$f_{\text{shock}}(U_s) = P_2(U_s) - P_{\text{eos}}\bigl(e_2(U_s),\, \rho_2(U_s)\bigr) = 0$$

in the interval $(v_{\text{piston}},\, U_{se})$.  Depending on the number and character of the roots, the solution falls into one of five regimes:

| # | Regime | Color | Description |
|---|--------|-------|-------------|
| 1 | **No plastic shock** | Blue | $v_p < v_Y$: piston velocity below elastic particle velocity; purely elastic response. |
| 2 | **Physical 2-wave** | Teal | Two roots; the weak-shock root (largest $U_s$, nearest $U_{se}$) gives moderate compression. Elastic precursor + plastic shock coexist. |
| 3 | **Strong shock only** | Yellow | Only the strong-shock root survives ($U_s / U_{se} \sim 0.15$--$0.3$, $\rho_2/\rho^Y \sim 1.9$--$2.0$). Unphysical for the piston problem. |
| 4 | **No root** | Red | $f_{\text{shock}}$ has no sign change; no self-consistent shock exists. |
| 5 | **Overdriven** | Purple | $v_p \geq U_{se}$: piston outruns the elastic precursor. |

## Regime maps — Mie-Gruneisen EOS

### Energy-split mode

![Energy-split regime map](../assets/regime_map_split.png)

### Original mode

![Original regime map](../assets/regime_map_orig.png)

### Comparison (both boundaries overlaid)

![Combined regime map](../assets/regime_map_combined.png)

## Regime maps — Ideal Gas EOS ($\gamma = 3$)

### Energy-split mode

![IG energy-split regime map](../assets/regime_map_ig_split.png)

### Original mode

![IG original regime map](../assets/regime_map_ig_orig.png)

### MG vs IG boundary comparison

![MG vs IG boundaries](../assets/regime_map_ig_vs_mg.png)

## $\log_{10}(U_{se}/U_s - 1)$ inside the physical 2-wave wedge

The quantity $U_{se}/U_s - 1$ measures the fractional speed excess of the elastic precursor over the plastic shock. Plotted on a $\log_{10}$ scale: $-3$ means the two waves differ by 0.1%, $-1$ by 10%, and $0$ by 100% (factor of 2). Grey areas are outside the physical 2-wave regime.

### 2x2 panel (all four cases)

![Ratio panel](../assets/ratio_panel.png)

### Individual maps

| | Energy-split | Original |
|---|---|---|
| **MG** | ![](../assets/ratio_mg_split.png) | ![](../assets/ratio_mg_orig.png) |
| **IG** | ![](../assets/ratio_ig_split.png) | ![](../assets/ratio_ig_orig.png) |

## Regime statistics

| Regime | MG split (%) | MG orig (%) | IG split (%) | IG orig (%) |
|--------|-------------|------------|-------------|-------------|
| No plastic shock | 34.4 | 35.1 | 22.3 | 23.7 |
| Physical 2-wave | 11.7 | 12.5 | 35.3 | 35.9 |
| Strong shock only | 46.2 | 45.2 | 0.0 | 0.0 |
| No root | 7.7 | 7.1 | 38.5 | 38.3 |
| Overdriven | 0.0 | 0.0 | 3.9 | 2.1 |

## Critical boundary: Physical 2-wave vs Strong shock only

The boundary between the physical 2-wave regime and the strong-shock-only regime is governed by the sign of $f_{\text{shock}}$ near $U_{se}$.

- When $f_{\text{shock}}(U_{se}^-) > 0$, the function starts positive near $U_{se}$ and crosses zero twice (weak + strong root).
- When $f_{\text{shock}}(U_{se}^-) < 0$, the weak-shock root has vanished and only the strong-shock root remains.

At the small-$Y_0/G$ limit, the critical piston velocity is approximately

$$v_p^* \approx 86879 \;\text{cm/s}\quad (v_p^*/C_0 \approx 0.163)$$

This threshold is an intrinsic property of the Mie-Gruneisen EOS ($C_0$, $s$, $\Gamma_0$), independent of strength.

The energy-split mode shifts this boundary to **lower** $v_{\text{piston}}$ because the EOS sees only $e_{\text{th}}$ instead of $e_{\text{total}}$. This lowers $P^Y$, which reduces the Hugoniot momentum pressure $P_2$ near $U_{se}$ more than it reduces $P_{\text{eos}}$, causing $f_{\text{shock}}$ to flip sign.

### Boundary table: maximum $v_p$ for physical 2-wave solution

| $Y_0/G$ | $v_p^*$ split (cm/s) | $v_p^*/C_0$ split | $v_p^*$ original (cm/s) | $v_p^*/C_0$ orig | Shift |
|---------|---------------------|----------------|------------------------|-----------------|-------|
| 0.005 | 86077 | 0.1615 | 86273 | 0.1619 | +0.23% |
| 0.015 | 85595 | 0.1606 | 86184 | 0.1617 | +0.68% |
| 0.025 | 85112 | 0.1597 | 86092 | 0.1615 | +1.14% |
| 0.035 | 84628 | 0.1588 | 86000 | 0.1614 | +1.59% |
| 0.045 | 84143 | 0.1579 | 85905 | 0.1612 | +2.05% |
| 0.055 | 83658 | 0.1570 | 85809 | 0.1610 | +2.51% |
| 0.065 | 83172 | 0.1560 | 85712 | 0.1608 | +2.96% |
| 0.075 | 82685 | 0.1551 | 85613 | 0.1606 | +3.42% |
| 0.085 | 82197 | 0.1542 | 85512 | 0.1604 | +3.88% |
| 0.095 | 81708 | 0.1533 | 85409 | 0.1602 | +4.33% |
| 0.105 | 81219 | 0.1524 | 85305 | 0.1600 | +4.79% |
| 0.115 | 80728 | 0.1515 | 85200 | 0.1598 | +5.25% |
| 0.125 | 80237 | 0.1505 | 85092 | 0.1596 | +5.71% |
| 0.135 | 79746 | 0.1496 | 84983 | 0.1594 | +6.16% |
| 0.145 | 79253 | 0.1487 | 84872 | 0.1592 | +6.62% |
| 0.155 | 78760 | 0.1478 | 84759 | 0.1590 | +7.08% |
| 0.165 | 78265 | 0.1468 | 84645 | 0.1588 | +7.54% |
| 0.175 | 77770 | 0.1459 | 84529 | 0.1586 | +8.00% |
| 0.185 | 77275 | 0.1450 | 84411 | 0.1584 | +8.45% |
| 0.195 | 76778 | 0.1440 | 84291 | 0.1581 | +8.91% |
| 0.205 | 76281 | 0.1431 | 84170 | 0.1579 | +9.37% |
| 0.215 | 75784 | 0.1422 | 84047 | 0.1577 | +9.83% |
| 0.225 | -- | -- | 83921 | 0.1575 | -- |
| 0.235 | -- | -- | 83794 | 0.1572 | -- |
| 0.245 | -- | -- | -- | -- | -- |
| 0.250 | -- | -- | -- | -- | -- |

## Key findings — Mie-Gruneisen

1. **The loss of the weak-shock root is tied to the Hugoniot reference curve** and its singularity at $\mu = 1/s$. Both modes lose the physical solution above $v_p / C_0 \approx 0.163$.

2. **Energy splitting shifts the critical boundary to lower $v_{\text{piston}}$.**  At $Y_0/G = 0.15$ the shift is large enough that $v_p = 80{,}000$ cm/s falls outside the physical 2-wave region in energy-split mode but remains inside it in the original mode.

3. **The physical 2-wave window is a triangular wedge** bounded below by $v_Y(Y_0/G)$ (diagonal line) and above by $v_p^*(Y_0/G)$ (near-horizontal curve). As $Y_0/G$ increases, $v_Y$ rises and $v_p^*$ drops, squeezing the wedge closed.

4. **At extreme velocities ($v_p / C_0 \gtrsim 0.65$), $f_{\text{shock}}$ loses all roots**, creating the red "no root" region. This is caused by the MG Hugoniot singularity ($\mu \to 1/s$).

## Key findings — Ideal Gas EOS

5. **The ideal gas EOS eliminates the Hugoniot singularity.** With $P = (\gamma - 1)\rho e$, there is no $(1 - s\mu)^2$ denominator. The pressure is well-defined for all densities.

6. **The "strong shock only" (yellow) region vanishes entirely.** Without the MG cold-curve contribution ($P_H$), the $f_{\text{shock}}$ function no longer develops the intermediate root structure that produces the unphysical strong-shock solution. The transition is directly from physical 2-wave to no-root.

7. **The "no root" region expands, not shrinks.** This is the key surprise. Without the cold pressure $P_H$, the EOS thermal pressure $(\gamma-1)\rho e$ exceeds the Rankine-Hugoniot momentum pressure $P_2$ across the entire $U_s$ range at moderate-to-high piston velocities. This is **not** a singularity artifact — the function is well-behaved (no NaN/inf) but simply never crosses zero. The two-wave structure genuinely cannot satisfy both the RH jump conditions and the ideal gas EOS simultaneously at those velocities.

8. **The physical 2-wave wedge is much larger** than with MG. The upper boundary ($v_p^*$) roughly doubles, extending to $v_p / C_0 \approx 0.3$ instead of $\approx 0.16$. This is because the ideal gas EOS lacks the stiff cold-curve pressure that flips $f_{\text{shock}}$ negative near $U_{se}$ in the MG case.

9. **The "Overdriven" (purple) region becomes visible** in the bottom-right corner, where $v_p$ exceeds $U_{se}$ at high $Y_0/G$. With MG the elastic precursor speed $U_{se}$ is much higher (due to $P_H$), pushing this region beyond the plotted range.

10. **The energy-split effect persists** with the ideal gas EOS, because elastic energy subtraction is EOS-independent. The split-mode 2-wave wedge is slightly smaller than the original-mode wedge, consistent with the MG behavior.

## Practical guidance

When using the energy-split formulation, ensure the piston velocity satisfies

$$v_Y(Y_0/G) < v_{\text{piston}} < v_p^*(Y_0/G)$$

The $v_p^*$ boundary can be pre-computed for a given material. The solver now issues a `warnings.warn()` when only the strong-shock root is found, alerting users that the result is unphysical.

The ideal gas EOS roughly doubles the physical 2-wave velocity window by removing the cold-curve pressure that causes the MG weak-shock root to vanish early. However, it replaces the MG "strong shock only" failure with a "no root" failure at similar velocities. Neither EOS produces valid two-wave solutions above $v_p / C_0 \approx 0.3$ (IG) or $\approx 0.16$ (MG). The ideal gas comparison confirms that the strong-shock-only artifact is specific to the MG cold curve, while the loss of *any* solution at high $v_p$ is a structural property of the two-wave elastoplastic model.
