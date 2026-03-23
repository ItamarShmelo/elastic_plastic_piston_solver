# Energy-Split Mode Comparison

This report compares the **original** (single total-energy) mode with the **energy-split** mode that decomposes internal energy into thermal and elastic parts: $e = e_{th} + e_{el}$, where $e_{el} = \boldsymbol{S}:\boldsymbol{S} \,/\, (4\rho G)$.

The EOS receives only the thermal energy $P(e_{th}, \rho)$ in energy-split mode. The difference grows with the ratio $Y_0/G$.

All cases use **CGS** units ($\mathrm{g/cm^3}$, $\mathrm{cm/s}$, $\mathrm{dyn/cm^2}$, $\mathrm{erg/g}$).

## Case 1: Aluminum

Material parameters: $\rho_0 = 2.79\,\mathrm{g/cm^3}$, $C_0 = 5.33 \times 10^5\,\mathrm{cm/s}$, $s = 1.34$, $\Gamma_0 = 2$, $G = 2.86 \times 10^{11}\,\mathrm{dyn/cm^2}$, $Y_0 = 2.6 \times 10^{9}\,\mathrm{dyn/cm^2}$, $v_{piston} = 5000\,\mathrm{cm/s}$.

$Y_0 / G = 0.0091$ — small elastic energy fraction.

| Quantity | Original | Energy-split | Rel. Diff (%) |
|----------|----------|--------------|---------------|
| $U_{se}$ | $652066$ | $651585$ | $0.0737$ |
| $U_s$ | $543900$ | $543892$ | $0.0013$ |
| $\rho^Y$ | $2.80271$ | $2.80271$ | $0.0000$ |
| $P^Y$ | $3.64661 \times 10^{9}$ | $3.63869 \times 10^{9}$ | $0.2174$ |
| $v^Y$ | $2957.21$ | $2955.03$ | $0.0737$ |
| $e^Y$ | $4.37255 \times 10^{6}$ | $4.3661 \times 10^{6}$ | $0.1473$ |
| $\rho_2$ | $2.81333$ | $2.81335$ | $0.0004$ |
| $P_2$ | $6.7437 \times 10^{9}$ | $6.73905 \times 10^{9}$ | $0.0690$ |
| $e_2$ | $1.37079 \times 10^{7}$ | $1.37031 \times 10^{7}$ | $0.0356$ |

Energy breakdown (split mode): $e_{th}^Y = 2.95841 \times 10^{6}\,\mathrm{erg/g}$, $e_{el}^Y = 1.4077 \times 10^{6}\,\mathrm{erg/g}$, $e_{th,2} = 1.22954 \times 10^{7}\,\mathrm{erg/g}$, $e_{el,2} = 1.4077 \times 10^{6}\,\mathrm{erg/g}$.

![Aluminum comparison](../assets/aluminum_comparison.png)

## Case 2: Copper

Material parameters: $\rho_0 = 8.93\,\mathrm{g/cm^3}$, $C_0 = 3.94 \times 10^5\,\mathrm{cm/s}$, $s = 1.49$, $\Gamma_0 = 2$, $G = 4.5 \times 10^{11}\,\mathrm{dyn/cm^2}$, $Y_0 = 9.0 \times 10^{8}\,\mathrm{dyn/cm^2}$, $v_{piston} = 2000\,\mathrm{cm/s}$.

$Y_0 / G = 0.002$ — very small elastic energy fraction.

| Quantity | Original | Energy-split | Rel. Diff (%) |
|----------|----------|--------------|---------------|
| $U_{se}$ | $472218$ | $472146$ | $0.0151$ |
| $U_s$ | $397696$ | $397695$ | $0.0001$ |
| $\rho^Y$ | $8.93893$ | $8.93893$ | $0.0000$ |
| $P^Y$ | $1.3903 \times 10^{9}$ | $1.3897 \times 10^{9}$ | $0.0432$ |
| $v^Y$ | $471.982$ | $471.91$ | $0.0151$ |
| $e^Y$ | $111383$ | $111350$ | $0.0302$ |
| $\rho_2$ | $8.97345$ | $8.97345$ | $0.0000$ |
| $P_2$ | $6.81592 \times 10^{9}$ | $6.81557 \times 10^{9}$ | $0.0051$ |
| $e_2$ | $2.1353 \times 10^{6}$ | $2.13516 \times 10^{6}$ | $0.0067$ |

Energy breakdown (split mode): $e_{th}^Y = 77777.5\,\mathrm{erg/g}$, $e_{el}^Y = 33572.2\,\mathrm{erg/g}$, $e_{th,2} = 2.10159 \times 10^{6}\,\mathrm{erg/g}$, $e_{el,2} = 33572.2\,\mathrm{erg/g}$.

![Copper comparison](../assets/copper_comparison.png)

## Case 3: High-Strength Material

To demonstrate a regime where the energy-split mode differs significantly, we use an artificial high-strength material with $Y_0 / G = 0.40$.

Material parameters: $\rho_0 = 2.79\,\mathrm{g/cm^3}$, $C_0 = 5.0 \times 10^4\,\mathrm{cm/s}$, $s = 1.2$, $\Gamma_0 = 2.5$, $G = 5.0 \times 10^{10}\,\mathrm{dyn/cm^2}$, $Y_0 = 2.0 \times 10^{10}\,\mathrm{dyn/cm^2}$, $v_{piston} = 4.0 \times 10^4\,\mathrm{cm/s}$.

| Quantity | Original | Energy-split | Rel. Diff (%) |
|----------|----------|--------------|---------------|
| $U_{se}$ | $201334$ | $175451$ | $12.8556$ |
| $U_s$ | $177574$ | $157943$ | $11.0551$ |
| $\rho^Y$ | $3.40771$ | $3.40771$ | $0.0000$ |
| $P^Y$ | $7.16708 \times 10^{9}$ | $2.23498 \times 10^{9}$ | $68.8160$ |
| $v^Y$ | $36495.7$ | $31803.9$ | $12.8556$ |
| $e^Y$ | $6.65967 \times 10^{8}$ | $5.05745 \times 10^{8}$ | $24.0585$ |
| $\rho_2$ | $3.49452$ | $3.64452$ | $4.2926$ |
| $P_2$ | $8.8518 \times 10^{9}$ | $5.75802 \times 10^{9}$ | $34.9509$ |
| $e_2$ | $8.2154 \times 10^{8}$ | $8.36182 \times 10^{8}$ | $1.7823$ |

Energy breakdown (split mode): $e_{th}^Y = 8.70329 \times 10^{7}\,\mathrm{erg/g}$, $e_{el}^Y = 4.18712 \times 10^{8}\,\mathrm{erg/g}$, $e_{th,2} = 4.1747 \times 10^{8}\,\mathrm{erg/g}$, $e_{el,2} = 4.18712 \times 10^{8}\,\mathrm{erg/g}$.

![High-strength comparison](../assets/high_strength_comparison.png)

## When Does the Energy-Split Mode Matter?

The elastic energy per unit mass at yield is $e_{el}^Y = Y_0^2 / (6 \rho^Y G)$. Its ratio to the total internal energy determines how much the two modes diverge. For typical metals ($Y_0/G \sim 10^{-3}$), the difference is negligible ($< 0.1\%$). The effect becomes significant ($> 1\%$) when $Y_0/G \gtrsim 0.05$, which can occur in very high-strength or ceramic-like materials.

