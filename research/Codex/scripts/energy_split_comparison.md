# Energy-Split Mode Comparison

This report compares the **original** (single total-energy) mode with the **energy-split** mode that decomposes internal energy into thermal and elastic parts: $e = e_{th} + e_{el}$, where $e_{el} = \boldsymbol{S}:\boldsymbol{S} \,/\, (4\rho G)$.

The EOS receives only the thermal energy $P(e_{th}, \rho)$ in energy-split mode. The difference grows with the ratio $Y_0/G$.

All cases use **CGS** units ($\mathrm{g/cm^3}$, $\mathrm{cm/s}$, $\mathrm{dyn/cm^2}$, $\mathrm{erg/g}$).

## Case 1: Aluminum

Material parameters: $\rho_0 = 2.79\,\mathrm{g/cm^3}$, $C_0 = 5.33 \times 10^5\,\mathrm{cm/s}$, $s = 1.34$, $\Gamma_0 = 2$, $G = 2.86 \times 10^{11}\,\mathrm{dyn/cm^2}$, $Y_0 = 2.6 \times 10^{9}\,\mathrm{dyn/cm^2}$, $v_{piston} = 5000\,\mathrm{cm/s}$.

$Y_0 / G = 0.0091$ — small elastic energy fraction.

| Quantity | Original | Energy-split | Rel. Diff (%) |
|----------|----------|--------------|---------------|
| $U_{se}$ | $652066$ | $651586$ | $0.0736$ |
| $U_s$ | $543900$ | $543895$ | $0.0008$ |
| $\rho^Y$ | $2.80271$ | $2.80271$ | $0.0000$ |
| $P^Y$ | $3.64661 \times 10^{9}$ | $3.6387 \times 10^{9}$ | $0.2170$ |
| $v^Y$ | $2957.21$ | $2955.03$ | $0.0736$ |
| $e^Y$ | $4.37255 \times 10^{6}$ | $4.36611 \times 10^{6}$ | $0.1471$ |
| $\rho_2$ | $2.81333$ | $2.81335$ | $0.0004$ |
| $P_2$ | $6.7437 \times 10^{9}$ | $6.73907 \times 10^{9}$ | $0.0686$ |
| $e_2$ | $1.37079 \times 10^{7}$ | $1.3703 \times 10^{7}$ | $0.0358$ |

Energy breakdown (split mode): $e_{th}^Y = 2.96055 \times 10^{6}\,\mathrm{erg/g}$, $e_{el}^Y = 1.40557 \times 10^{6}\,\mathrm{erg/g}$, $e_{th,2} = 1.23028 \times 10^{7}\,\mathrm{erg/g}$, $e_{el,2} = 1.40025 \times 10^{6}\,\mathrm{erg/g}$.

![Aluminum comparison](../assets/aluminum_comparison.png)

## Case 2: Copper

Material parameters: $\rho_0 = 8.93\,\mathrm{g/cm^3}$, $C_0 = 3.94 \times 10^5\,\mathrm{cm/s}$, $s = 1.49$, $\Gamma_0 = 2$, $G = 4.5 \times 10^{11}\,\mathrm{dyn/cm^2}$, $Y_0 = 9.0 \times 10^{8}\,\mathrm{dyn/cm^2}$, $v_{piston} = 2000\,\mathrm{cm/s}$.

$Y_0 / G = 0.002$ — very small elastic energy fraction.

| Quantity | Original | Energy-split | Rel. Diff (%) |
|----------|----------|--------------|---------------|
| $U_{se}$ | $472218$ | $472146$ | $0.0151$ |
| $U_s$ | $397696$ | $397696$ | $0.0000$ |
| $\rho^Y$ | $8.93893$ | $8.93893$ | $0.0000$ |
| $P^Y$ | $1.3903 \times 10^{9}$ | $1.3897 \times 10^{9}$ | $0.0432$ |
| $v^Y$ | $471.982$ | $471.91$ | $0.0151$ |
| $e^Y$ | $111383$ | $111350$ | $0.0302$ |
| $\rho_2$ | $8.97345$ | $8.97345$ | $0.0000$ |
| $P_2$ | $6.81592 \times 10^{9}$ | $6.81557 \times 10^{9}$ | $0.0051$ |
| $e_2$ | $2.1353 \times 10^{6}$ | $2.13516 \times 10^{6}$ | $0.0067$ |

Energy breakdown (split mode): $e_{th}^Y = 77788.7\,\mathrm{erg/g}$, $e_{el}^Y = 33561\,\mathrm{erg/g}$, $e_{th,2} = 2.10173 \times 10^{6}\,\mathrm{erg/g}$, $e_{el,2} = 33431.9\,\mathrm{erg/g}$.

![Copper comparison](../assets/copper_comparison.png)

## Case 3: High-Strength Material

To demonstrate a regime where the energy-split mode differs significantly, we use an artificial high-strength material with $Y_0 / G = 0.70$.

Material parameters: $\rho_0 = 2.79\,\mathrm{g/cm^3}$, $C_0 = 5.33 \times 10^5\,\mathrm{cm/s}$, $s = 1.34$, $\Gamma_0 = 2$, $G = 2.86 \times 10^{11}\,\mathrm{dyn/cm^2}$, $Y_0 = 2.0 \times 10^{11}\,\mathrm{dyn/cm^2}$, $v_{piston} = 4.5 \times 10^5\,\mathrm{cm/s}$.

| Quantity | Original | Energy-split | Rel. Diff (%) |
|----------|----------|--------------|---------------|
| $U_{se}$ | $1.02747 \times 10^{6}$ | $978930$ | $4.7246$ |
| $U_s$ | $752030$ | $797190$ | $6.0050$ |
| $\rho^Y$ | $3.95781$ | $3.95781$ | $0.0000$ |
| $P^Y$ | $7.35756 \times 10^{11}$ | $6.55574 \times 10^{11}$ | $10.8979$ |
| $v^Y$ | $303172$ | $288848$ | $4.7246$ |
| $e^Y$ | $4.59567 \times 10^{10}$ | $4.17167 \times 10^{10}$ | $9.2260$ |
| $\rho_2$ | $5.88185$ | $5.79487$ | $1.4788$ |
| $P_2$ | $9.96595 \times 10^{11}$ | $9.79798 \times 10^{11}$ | $1.6854$ |
| $e_2$ | $1.28566 \times 10^{11}$ | $1.17892 \times 10^{11}$ | $8.3027$ |

Energy breakdown (split mode): $e_{th}^Y = 3.58271 \times 10^{10}\,\mathrm{erg/g}$, $e_{el}^Y = 5.88962 \times 10^{9}\,\mathrm{erg/g}$, $e_{th,2} = 1.13869 \times 10^{11}\,\mathrm{erg/g}$, $e_{el,2} = 4.02252 \times 10^{9}\,\mathrm{erg/g}$.

![High-strength comparison](../assets/high_strength_comparison.png)

## When Does the Energy-Split Mode Matter?

The elastic energy per unit mass at yield is $e_{el}^Y = Y_0^2 / (6 \rho^Y G)$. Its ratio to the total internal energy determines how much the two modes diverge. For typical metals ($Y_0/G \sim 10^{-3}$), the difference is negligible ($< 0.1\%$). The effect becomes significant ($> 1\%$) when $Y_0/G \gtrsim 0.05$, which can occur in very high-strength or ceramic-like materials.

