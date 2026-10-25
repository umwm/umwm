# UMWM Technical Reference

This file is a technical reference for the physics implemented in the University of Miami Wave Model (UMWM).

## Implementation Map

| Physics area | Main implementation | Purpose |
| --- | --- | --- |
| Model domain and loop | `src/umwm_domain.f90` | `domain_type` owns the model state, initializes it, and runs the time loop: source functions, source integration, advection, refraction, stress, Stokes drift, diagnostics, and output. |
| Output and output-time physics | `src/umwm_domain_output.f90` | Integrated diagnostics, ocean-side stress and energy fluxes, gridded and spectrum NetCDF output. |
| Forcing state | `src/umwm_forcing.f90` | `forcing_type` owns mutable wind, current, density, sea-ice, remapped, and interpolation-scratch fields. |
| Source functions | `src/umwm_source_functions.f90` | Wind input, wave breaking, nonlinear downshifting, turbulence, and sea-ice attenuation. |
| Source integration | `src/umwm_physics.f90` | Exponential source update, diagnostic spectral tail, physics time step limit. |
| Propagation and refraction | `src/umwm_advection.f90` | First-order upstream advection in geographic and directional space. |
| Stress | `src/umwm_stress.f90` | Form/skin drag, drag coefficient, stress tail, and spectral stress integrals. The domain evaluates the atmospheric stress each time step, and the ocean-top and ocean-bottom momentum fluxes (`ocean_stress` in `src/umwm_domain_output.f90`) at gridded output times, or each time step when coupled with ESMF. |
| Stokes drift | `src/umwm_stokes.f90` | `stokes_type` with wave-induced Stokes drift and e-folding depth. |
| Spectral grid | `src/umwm_spectrum.f90` | `spectrum_type`, frequency and direction coordinates, and spectral spacings. |
| Dispersion | `src/umwm_dispersion.f90` | Capillary-gravity dispersion relation, wavenumber solve, angular frequency, and group speed. |
| Initialization and precomputed factors | `src/umwm_domain.f90` | Runtime allocation, spectral coordinates, wave kinematics, integration weights, dissipation constants, CFL limits. |
| Parallel environment and halo exchange | `src/umwm_env.f90`, `src/umwm_mpi.f90` | Process rank and count; halo exchange and gathering of tile fields. |
| Command line | `src/umwm_cli.f90` | Arguments of the `umwm` executable: the main namelist to run, by path or as an example case of `examples/` (`--case`), and the check of the run directory. |
| C API | `src/umwm_capi.f90` | `bind(c)` procedures for foreign-language callers: the domain lifecycle through an opaque handle, and access to the domain fields by name as pointers into the model memory. Used by the Python package in `python/`. |

The main program instantiates one `domain_type` from `src/umwm_domain.f90`
and calls its `initialize`, `run`, and `finalize` procedures, with the
main namelist given on the command line (`namelists/main.nml` by default);
the `spectrum.nml` and `exclude.nml` namelists are read from the directory
of the main namelist. `run` loops over `advance`, which integrates one
global time step `dtg`: the forcing update, the inner loop of source time
steps, and the output due at the new time; `advance` and `finished` are
public so that a driver, such as the Python package, can step the model
and inspect it between steps. The
domain is constructed from `config_type`, `spectrum_type`, `grid_type`,
`forcing_type`, and `stokes_type` components, and owns the wave spectrum
and all fields derived from it; there is no global model state.

The Python package in `python/` (see `python/README.md`) drives the model
through the C API of `src/umwm_capi.f90`, compiled into the shared library
`libumwm.so` by `make shared`: a Python `Domain` holds a `domain_type`
through an opaque handle, and the fields of the domain are exposed as
NumPy arrays viewing the Fortran memory, so the time integration runs
entirely in Fortran.
`config_type` remains the source of forcing switches and constants, while
`forcing_type` owns the mutable runtime forcing snapshots and remapped
fields.

The source functions, integration, advection, and stress kernels operate
on one spectral column (all frequencies and directions of a grid cell)
with explicit-shape array arguments, so that compilers can vectorize them;
the domain loops over cells and passes column slices of its state.

Within each global forcing/output time step `dtg`, UMWM repeatedly runs:

1. Advance file-backed forcing snapshots, when not ESMF-coupled.
2. Interpolate `forcing` fields.
3. Compute `Sin`, `Sds`, `Sice`, turbulence and linear sink rates, and the
   largest stable physics time step of each cell.
4. Find the dynamic physics time step `dts`, the minimum over all cells.
5. Integrate source terms over `dts`, including `Snl`.
6. Exchange MPI halos when enabled.
7. Propagate wave variance geographically.
8. Refract wave variance directionally.
9. Compute atmospheric stress and drag.
10. At output times, compute Stokes drift, diagnostics, ocean stress, spectra, grid output, and restart output.
    When coupled with ESMF, ocean stress is computed each time step.

Steps 5 and 7 to 9 are pipelined over grid cells: a cell is propagated a
fixed number of cells after its integration, at least the largest
distance to a neighbor cell in the tile, so that the spectra of its
neighbors are integrated and still in cache. Under MPI, the cells in the
halo of neighboring tiles are integrated first, and the halo is exchanged
before the pipeline starts. The halo is one grid row or column wide, and
the model stops at startup if a tile is narrower than the halo of its
neighbors.

## State Variables and Conventions

UMWM predicts the surface elevation variance spectrum
$E(\mathbf{x}, k, \phi, t)$ on a two-dimensional horizontal grid. The
state is discretized in:

- horizontal grid cell index $i$,
- logarithmically spaced frequency bins $f_o$,
- corresponding wavenumber bins $k_o(i)$,
- uniformly spaced propagation direction bins $\phi_p$,
- time.

The main spectrum array is `e(o,p,i)`.

Forcing fields are represented by `type(forcing_type)` from
`src/umwm_forcing.f90`. It stores native-grid snapshots for backward and
forward forcing time levels, interpolated native-grid scratch fields,
and remapped one-dimensional fields consumed by source functions,
advection/refraction, stress, diagnostics, and output. File-backed
updates remain disabled under `ESMF`, where coupled forcing is expected
to be supplied externally.

The spectral grid metadata is represented explicitly by
`type(spectrum_type)` from `src/umwm_spectrum.f90`. `domain % initialize`
constructs one spectrum object after namelist input is read, and the
domain uses it for the spectral dimensions and spacings, the direction
projection of the grid, restart I/O, and spectrum and grid output. The
domain keeps default-real copies of the frequency and direction
coordinates, `f(:)` and `th(:)`, with their cosines and sines, for the
model arrays and formulas.

The model uses mathematical direction convention: $\phi=0$ points in the
positive $x$ direction and positive angles rotate counter-clockwise.
Wind direction `wdir` uses the same convention.

The relationship between wave energy spectrum $E'$ and variance spectrum
$E$ is

$$
E'(\mathbf{x}, k, \phi) = \rho_w g E(\mathbf{x}, k, \phi).
$$

For implementation, the variance spectrum is the prognostic quantity.
Most integrated quantities use the discrete spectral weight

$$
k \, dk \, d\phi.
$$

In code, `dth` is $d\phi$, `dwn(o,i)` is $dk$, and `kdk(o,i)` is
$k_o \, dk_o$.

## Spectral Grid and Wave Kinematics

Frequencies are logarithmically spaced:

$$
\Delta \ln f =
\frac{\ln f_{\max} - \ln f_{\min}}{N_f - 1},
\qquad
f_o = \exp\left(\ln f_{\min} + (o-1)\Delta\ln f\right).
$$

This construction is implemented by `spectrum_type` using the `DOMAIN`
namelist values `om`, `pm`, `fmin`, and `fmax`:

| `spectrum_type` field | Meaning |
| --- | --- |
| `num_frequencies` | Number of frequency bins, initialized from `om`. |
| `num_directions` | Number of direction bins, initialized from `pm`. |
| `frequency_min`, `frequency_max` | Lowest and highest frequency bins, initialized from `fmin` and `fmax`. |
| `frequency(:)` | Logarithmically spaced frequency coordinates. |
| `direction(:)` | Uniformly spaced directional coordinates. |
| `dlnf` | Log-frequency spacing. |
| `dth` | Directional spacing. |

Directions are uniformly spaced over a full circle:

$$
\Delta\phi = \frac{2\pi}{N_\phi},
\qquad
\phi_p =
\left(p-\frac{N_\phi+1}{2}\right)\Delta\phi.
$$

During initialization, the domain arrays `f(:)` and `th(:)` are assigned
from `spectrum % frequency` and `spectrum % direction`.

The angular-frequency increment used to infer $dk$ is

$$
d\omega_o = 2\pi \Delta \ln f \, f_o.
$$

For every grid cell and frequency, the wavenumber is solved from the
finite-depth capillary-gravity dispersion relation

$$
\omega^2 =
\left(gk + \frac{\sigma}{\rho_w} k^3\right)\tanh(kd),
\qquad
\omega = 2\pi f,
$$

where $d$ is water depth and $\sigma$ is water surface tension
(`sfct`). The pure elemental function `wavenumber()` in
`src/umwm_dispersion.f90` solves a nondimensional form by Newton
iteration. The same module also provides `angular_frequency()` for the
forward relation and `group_speed()` for the derivative used by the
model.

Intrinsic phase speed and group speed are

$$
c = \frac{\omega}{k},
$$

$$
c_g =
c\left[
\frac{1}{2}
+ \frac{\min(kd,20)}{\sinh(2\min(kd,20))}
+ \frac{\sigma k^2}{\rho_w g + \sigma k^2}
\right].
$$

The `min(kd,20)` cap is the implementation limit used in
`group_speed()` to avoid hyperbolic overflow at large nondimensional
depth.

The spectral wavenumber increment is approximated from the group speed:

$$
dk_o = \frac{d\omega_o}{|c_{g,o}|}.
$$

Precomputed factors include:

$$
kdk = k\,dk,\qquad
k3dk = k^3\,dk,\qquad
k4 = k^4,\qquad
fkovg = \frac{fk}{g}.
$$

The half wavelength used by the wind input is

$$
\frac{\lambda}{2} = \frac{|c|}{2f}.
$$

The logarithmic wind-profile factor is capped at the field-log-layer top:

$$
\log\left(\frac{\lambda/2}{z}\right)
\rightarrow
\log\left(\frac{20\,\mathrm{m}}{z}\right)
\quad \text{when } \lambda/2 > 20\,\mathrm{m}.
$$

## Governing Equation

The model solves the wave energy balance equation:

$$
\frac{\partial E'}{\partial t}
+ \frac{\partial(\dot{\mathbf{x}}E')}{\partial \mathbf{x}}
+ \frac{\partial(\dot{k}E')}{\partial k}
+ \frac{\partial(\dot{\phi}E')}{\partial \phi}
=
\rho_w g \sum_{i=1}^{n} S_i.
$$

The geographic propagation velocity is

$$
\dot{\mathbf{x}} = \mathbf{c}_g + \mathbf{u},
$$

where $\mathbf{u}=(u,v)$ is the current in the wave boundary layer.

Wavenumber-space advection $\dot{k}$ is nonzero for changing currents
or moving bottoms. The implemented model treats currents as
quasi-stationary on wave growth/decay time scales and neglects this
term.

Directional-space advection is refraction. It is nonzero for spatially
varying depth and/or inhomogeneous currents.

For the variance spectrum, the local source terms are written as signed
rates or tendencies acting on $E$.

## Source-Term Sign Conventions

The manual writes physical source terms as signed variance-spectrum
tendencies. The implementation stores several of them as positive decay
rates. This distinction is important.

| Physical term | Manual sign | Implementation array | Implementation meaning |
| --- | --- | --- | --- |
| Wind input/export | positive or negative | `ssin(o,p,i)` | Signed rate multiplying `e`. |
| Breaking dissipation | negative | `sds(o,p,i)` | Positive decay rate multiplying `e`. |
| Turbulent dissipation | negative | `sdt(o,i)` | Positive decay rate multiplying `e`. |
| Viscous dissipation | negative | `sdv(o,i)` | Positive decay rate multiplying `e`. |
| Bottom friction and percolation | negative | `sbf(o,i)` | Positive combined bottom decay rate multiplying `e`. |
| Nonlinear downshifting | source/sink | `snl(o,p,i)` | Absolute tendency, already multiplied by donor-bin variance. |
| Sea-ice attenuation | negative when active | `sice(o,i)` | Signed rate, normally negative. |

For prognostic bins, the source update uses

$$
E_s^{n+1}
=
E^n
\exp\left[
\Delta t
\left(
S_{in}^*
- S_{ds}^*
- S_{bottom}^*
- S_{dt}^*
- S_{dv}^*
+ S_{ice}^*
\right)
\right]
+ \Delta t\,S_{nl},
$$

where stars denote the implementation rate arrays that multiply $E$.
`Snl` is not inside the exponential because it is carried as an
absolute bin-to-bin tendency.

## Wind Input `Sin`

Wind input follows Jeffreys' sheltering hypothesis. Define the wave-wind
angle

$$
\theta = \psi - \phi,
$$

where $\psi$ is wind direction, and define the effective wind-wave
relative speed

$$
\Delta U =
U_{\lambda/2}\cos\theta
- c
- u\cos\phi
- v\sin\phi.
$$

The physical source is

$$
S_{in}
=
A_1
\Delta U |\Delta U|
\frac{k\omega}{g}
\frac{\rho_a}{\rho_w}
E(k,\phi).
$$

The corresponding implementation rate is

$$
S_{in}^*
=
A_1
\Delta U |\Delta U|
\frac{k\omega}{g}
\frac{\rho_a}{\rho_w}.
$$

The wind speed entering this term is evaluated at half wavelength above
the surface:

$$
U_{\lambda/2}
=
U_z
+ \frac{u_*}{\kappa}
\log\left(\frac{\lambda/2}{z}\right),
$$

with $\kappa=0.4$, hence $1/\kappa=2.5$ in the code. In coupled `ESMF`
builds, stability-function corrections are added to this log-profile
term.

The sheltering coefficient depends on sea state:

$$
A_1 =
\begin{cases}
\text{variable}, & \Delta U > 0 \quad \text{wind sea},\\
0.001, & 0 < U_{\lambda/2}\cos\theta < c+u\cos\phi+v\sin\phi \quad \text{swell with wind},\\
0.1, & \cos\theta < 0 \quad \text{swell against wind}.
\end{cases}
$$

In the current source, the variable positive-input coefficient is
computed by `sheltering_coare35(wspd)`, then negative input is scaled by
`sin_diss1/sin_fac`, and swell-overrunning-wind input is further scaled
by `sin_diss2/sin_diss1`.

The `sheltering_coare35` approximation is piecewise in wind speed
$U$:

$$
A_1(U) =
\begin{cases}
i_0 + m_1 U, & U \le x_1,\\
a + bU + cU^2, & x_1 < U \le x_2,\\
s_1\exp\left[-\dfrac{U-x_2}{d_c U}\right], & U > x_2,
\end{cases}
$$

with

$$
\begin{aligned}
x_1 &= 15,\quad x_2 = 33,\quad x_3 = 60,\\
y_1 &= 0.10,\quad y_2 = 0.09,\quad y_3 = 0.06,\\
i_0 &= 0.04,\quad c_r = 0.65,\quad d_c = 1.6,\\
m_1 &= \frac{y_1-i_0}{x_1},\quad
m_2 = \frac{y_3-y_2}{x_3-x_2},\\
c &= c_r\frac{m_2-m_1}{x_2-x_1},\\
b &= m_1 - 2cx_1,\\
a &= y_1 - bx_1 - cx_1^2,\\
s_1 &= a + bx_2 + cx_2^2.
\end{aligned}
$$

Sea ice reduces wind input by multiplying `ssin` by $(1-f_{ice})$.
For bins above the local diagnostic cutoff, negative wind input is
clipped to zero.

## Breaking Dissipation `Sds`

Breaking dissipation is nonlinear in the saturation spectrum

$$
B(k,\phi) = k^4 E(k,\phi).
$$

The physical dissipation source is

$$
S_{ds}(k,\phi)
=
-A_2
\coth(0.2kd)
\left(1 + A_3\overline{\chi^2}(k,\phi)\right)^2
B^n(k,\phi)\,
\omega E(k,\phi).
$$

The implementation stores the positive decay rate

$$
S_{ds}^*
=
A_2
\coth(0.2kd)
\left(1 + A_3\overline{\chi^2}\right)^2
B^n
\omega.
$$

Default coefficients are

$$
A_2 = 42,\qquad
A_3 = 360,\qquad
n = 2.4.
$$

The longer-wave mean-square-slope contribution is computed directionally:

$$
\overline{\chi^2}_{o,p}
=
\sum_{q<o}\sum_{p'}
E_{q,p'}
\cos^2(\phi_{p'}-\phi_p)
k_q^3\,dk_q\,d\phi.
$$

Because $\cos^2 x = (1 + \cos 2x)/2$, the double sum reduces to three
direction-integrated moments of the longer waves,

$$
\overline{\chi^2}_{o,p}
=
M^{(0)}_o + \cos(2\phi_p)\,M^{(c)}_o + \sin(2\phi_p)\,M^{(s)}_o,
\qquad
M^{(0,c,s)}_o
=
\frac{d\phi}{2}
\sum_{q<o}
k_q^3\,dk_q
\sum_{p'}
E_{q,p'}\,\{1, \cos 2\phi_{p'}, \sin 2\phi_{p'}\},
$$

which `mss_moments` evaluates in $O(N_f N_\phi)$ rather than
$O(N_f N_\phi^2)$ operations. The power of the saturation spectrum is
evaluated in double precision.

The code computes and stores the spilling part, before applying
$\coth(0.2kd)$. This lets `Snl` use the spilling-only dissipation
rate. The factor `cothkd` is applied where `sds` is used, to include the
finite-depth/plunging-breaker enhancement in the decay rate used by the
source integrator, the stress diagnostics, and the `Sds` output. Between
output times, `sds` is evaluated only for the frequencies used by the
integration, up to two bins above the prognostic range.

## Turbulent Dissipation `Sdt`

Ambient turbulence in the wave boundary layer attenuates waves:

$$
S_{dt}
=
-A_4 u_{*w}kE(k,\phi).
$$

The implementation stores the positive decay rate

$$
S_{dt}^*
=
A_4 u_{*w} k.
$$

The code approximates the water-side friction velocity from the air-side
friction velocity and density ratio:

$$
u_{*w}
=
u_*\sqrt{\frac{\rho_a}{\rho_w}}.
$$

Default:

$$
A_4 = 0.002.
$$

## Viscous Dissipation `Sdv`

Viscous dissipation converts wave energy to heat:

$$
S_{dv}
=
-4\nu k^2 E(k,\phi).
$$

The implementation stores

$$
S_{dv}^*
=
4\nu k^2,
$$

where $\nu$ is `nu_water`. This term is usually important only for the
shortest waves in clean water.

## Bottom Friction and Percolation `Sbf`

The manual separates bottom friction and bottom percolation:

$$
S_{bf}
=
-G_f\frac{k}{\sinh(2kd)}E(k,\phi),
$$

$$
S_{bp}
=
-G_p\frac{k}{\cosh^2(kd)}E(k,\phi).
$$

The implementation combines both into the positive decay-rate array
`sbf`:

$$
S_{bottom}^*
=
G_f\frac{k}{\sinh(2kd)}
+
G_p\frac{k}{\cosh^2(kd)}.
$$

Typical coefficient ranges are:

$$
G_f = 0.001\text{ to }0.01\,\mathrm{m\,s^{-1}},
\qquad
G_p = 0.0006\text{ to }0.01\,\mathrm{m\,s^{-1}}.
$$

Defaults in `namelists/main.nml` are `sbf_fac = 0.003` and
`sbp_fac = 0.003`.

## Nonlinear Downshifting `Snl`

UMWM represents nonlinear wave-wave transfer by moving a quantity
proportional to spilling dissipation into the next two lower-wavenumber
bins. The manual form is

$$
S_{nl}(k,\phi)
=
A_5
\left[
b_1 S_{sb}(k-\Delta k,\phi)
+ b_2 S_{sb}(k-2\Delta k,\phi)
- S_{sb}(k,\phi)
\right],
$$

where

$$
b_1 = \exp\left[-16\left(\frac{\Delta f}{f}\right)^2\right],
\qquad
b_2 = \exp\left[-16\left(2\frac{\Delta f}{f}\right)^2\right],
$$

renormalized so $b_1+b_2=1$, and

$$
S_{sb} = \frac{S_{ds}}{\coth(0.2kd)}
$$

is the spilling-only breaking term.

In the discrete implementation, bin $o$ receives transfer from bins
$o+1$ and $o+2$ and loses energy to lower-wavenumber bins:

$$
S_{nl,o,p}
=
\beta_{1,o} S_{ds,o+1,p}^*E_{o+1,p}
+
\beta_{2,o} S_{ds,o+2,p}^*E_{o+2,p}
-
A_5 S_{ds,o,p}^*E_{o,p}.
$$

The implementation factors are

$$
\beta_{1,o}
=
A_5 b_1
\frac{(k\,dk)_{o+1}}{(k\,dk)_o},
\qquad
\beta_{2,o}
=
A_5 b_2
\frac{(k\,dk)_{o+2}}{(k\,dk)_o}.
$$

`snl_fac` is the implementation coefficient $A_5$; the default is `5.0`.
The array `snl_arg = 1 - (\beta_1+\beta_2)` is used in the dynamic
physics time-step estimate to account for the nonlinear transfer's
effective sink strength.

## Sea-Ice Attenuation `Sice`

The manual identifies sea-ice attenuation as an active source/sink term.
The current source implementation follows the attenuation form cited in
the source comments as Kohout et al. (2014).

If the sea-ice fraction is below `fice_lth`, `sice = 0`. If it is above
`fice_lth`, the code computes local significant wave height

$$
H_s = 4\sqrt{\sum_{o,p} E_{o,p} k_o\,dk_o\,d\phi}.
$$

Then it defines an attenuation coefficient

$$
\alpha =
\begin{cases}
C_1 H_s, & H_s < H_{th},\\
C_1 H_{th}, & H_s \ge H_{th},
\end{cases}
$$

with

$$
H_{th}=3\,\mathrm{m},
\qquad
C_1=-5.35\times10^{-6}\,\mathrm{m^{-1}}.
$$

The implemented rate is

$$
S_{ice}^* = 2c_g\alpha.
$$

Because $C_1<0$, this is normally a sink. When `fice > fice_uth`,
propagation and refraction set the updated spectrum to zero in that
cell. The default thresholds are `fice_lth = 0.30` and
`fice_uth = 0.75`.

## Source Integration and Diagnostic Tail

The local source contribution is split from propagation and refraction:

$$
\frac{\partial E}{\partial t}
=
\left(\frac{\partial E}{\partial t}\right)_s
+
\left(\frac{\partial E}{\partial t}\right)_a.
$$

For source rates that multiply $E$,

$$
\left(\frac{\partial E}{\partial t}\right)_s
=
\sum_i S_i^*E.
$$

The exponential source solution is

$$
E_s^{n+1}
=
E^n
\exp\left(\sum_i S_i^*\Delta t\right).
$$

The timestep is limited so the exponential factor remains bounded:

$$
\exp\left(\sum_i S_i^*\Delta t\right) < r.
$$

In code, `explim` is the logarithmic exponent limit and the diagnostic
physics timestep is

$$
\Delta t_{phys}(i)
=
\frac{\texttt{explim}}
{\max_{o,p}|\Lambda_{o,p,i}|},
$$

where the estimated local rate is

$$
\Lambda
=
S_{in}^*
- S_{ds}^*
- S_{bottom}^*
- S_{dt}^*
- S_{dv}^*
+ S_{ice}^*.
$$

The actual source timestep is

$$
dts =
\min\left(
\min_i \Delta t_{phys}(i),
dtamin,
dtg-sumt
\right).
$$

`dtamin` is the advective CFL limit computed at initialization:

$$
dtamin = 0.98\,cfl_{lim}
\frac{\min(\Delta x,\Delta y)}{\max(c_g)}.
$$

The CFL factor is

$$
cfl_{lim} =
\begin{cases}
\cos\left(\frac{\pi}{4}-\frac{\Delta\phi}{2}\right), & N_\phi \text{ divisible by } 8,\\
\frac{1}{\sqrt{2}}, & \text{otherwise}.
\end{cases}
$$

For the first non-restart step, `dts` is set to zero so only the
diagnostic spectrum is initialized.

The prognostic/diagnostic cutoff frequency is

$$
f_c = \texttt{fcut\_fac}\,\frac{g}{U_{10}},
$$

by default with `fcut_fac = 0.53`, so that $f_c \approx 4f_{PM}$, and
limited by `fprog`. Where waves are far younger than a fully developed
sea, as in a laboratory tank, their spectral peak lies above this default
cutoff, and a larger `fcut_fac` keeps them in the prognostic range. The
local cutoff bin `oc(i)` marks the highest prognostic bin.

For bins above `oc(i)`, UMWM assumes a wind-aligned equilibrium range.
When the numerator is nonnegative, the diagnostic spectrum is assigned
from the balance of wind input, turbulence, viscosity, sea ice, and
breaking:

$$
E_{eq}
=
k^{-4}
\left[
\frac{
S_{in}^* - S_{dt}^* - S_{dv}^* + S_{ice}^*
}{
2\pi A_2 f
\left(1+A_3\overline{\chi^2}\right)^2
\coth(0.2kd)
}
\right]^{1/n}.
$$

After the source calculation, the implementation uses the time-split
intermediate state

$$
E^* = \frac{E^n + E_s^{n+1}}{2}
$$

for propagation and refraction.

## Geographic Propagation

The advection equation in Cartesian projection is

$$
\frac{\partial E}{\partial t}
=
-\frac{\partial[(c_g\cos\phi+u)E]}{\partial x}
-\frac{\partial[(c_g\sin\phi+v)E]}{\partial y}
-\frac{\partial(\dot{\phi}E)}{\partial \phi}.
$$

Geographic propagation uses first-order upstream differencing. For one
dimension,

$$
\frac{\partial(\dot{x}E)}{\partial x}
\approx
\frac{\Phi_{i+1/2}-\Phi_{i-1/2}}{\Delta x},
$$

with

$$
\Phi_{i+1/2}
=
\frac{\dot{x}_{i+1/2}+|\dot{x}_{i+1/2}|}{2}E_i
+
\frac{\dot{x}_{i+1/2}-|\dot{x}_{i+1/2}|}{2}E_{i+1},
$$

$$
\Phi_{i-1/2}
=
\frac{\dot{x}_{i-1/2}+|\dot{x}_{i-1/2}|}{2}E_{i-1}
+
\frac{\dot{x}_{i-1/2}-|\dot{x}_{i-1/2}|}{2}E_i,
$$

and

$$
\dot{x}_{i\pm1/2}
=
\frac{\dot{x}_i+\dot{x}_{i\pm1}}{2}.
$$

The implementation applies this form over east, west, north, and south
faces using cell-face lengths and cell area. Propagation by currents is
added only when the current fields are nonzero.

Open boundary conditions are applied next to land and at regional domain
edges unless boundary input is supplied. Global domains use periodic
east-west boundary conditions. Land cells take, at each frequency, the
smallest phase and group speed over all sea cells.

## Refraction

The directional rotation rate is

$$
\dot{\phi}
=
\frac{\partial(c\sin\phi+v)}{\partial x}
-
\frac{\partial(c\cos\phi+u)}{\partial y}.
$$

Positive $\dot{\phi}$ rotates wave energy counter-clockwise; negative
$\dot{\phi}$ rotates it clockwise.

Refraction is discretized with the same upstream-flux form, replacing
$x$ by $\phi$. Directional space is periodic, so no directional boundary
condition is required.

The stability condition is

$$
\frac{\dot{\phi}\Delta t}{\Delta\phi} < 1.
$$

The implementation computes a refraction timestep `dtr` no larger than
`dts`; if needed, the directional rotation is effectively limited to the
stable step.

## Wind Stress and Drag

Wave energy and momentum are related by phase speed:

$$
E = Mc.
$$

The resolved wind form stress is computed from the wind-input source:

$$
\tau_x
=
\rho_w g
\int_{-\pi}^{\pi}
\int_{k_{\min}}^{k_{\max}}
\frac{S_{in}}{c}
\cos\phi\,k\,dk\,d\phi,
$$

$$
\tau_y
=
\rho_w g
\int_{-\pi}^{\pi}
\int_{k_{\min}}^{k_{\max}}
\frac{S_{in}}{c}
\sin\phi\,k\,dk\,d\phi.
$$

The form drag coefficient is

$$
C_{d,f}
=
\frac{\sqrt{\tau_x^2+\tau_y^2}}{\rho_a U_z^2}.
$$

The smooth-flow skin drag coefficient is computed from

$$
C_{d,s}
=
\frac{u_*^2}{U^2(z)}
=
\frac{\kappa^2}{\left[\log(z/z_0)\right]^2},
$$

with molecular-sub-layer roughness

$$
z_0 = 0.132\frac{\nu_a}{u_*}.
$$

The implementation solves this relation by six fixed-point iterations.
It uses wind relative to the moving surface:

$$
\mathbf{U}_{rel}
=
\mathbf{U}_{wind}
- \left(\mathbf{u}_{current}+\mathbf{u}_{Stokes}(z=0)\right).
$$

The skin drag is reduced by wave sheltering:

$$
C_{d,s}^{new}
=
\frac{C_{d,s}^{old}}{3}
\left(
1 + 2\frac{C_{d,s}^{old}}{C_{d,s}^{old}+C_{d,f}}
\right),
$$

and capped at $10^{-2}$.

The total atmospheric stress and drag coefficient are

$$
\boldsymbol{\tau}
=
\boldsymbol{\tau}_{form}
+
\boldsymbol{\tau}_{skin},
\qquad
C_d
=
\frac{|\boldsymbol{\tau}|}{\rho_a U_z^2}.
$$

The air-side friction velocity is updated from the total stress:

$$
u_*
=
\sqrt{\frac{|\boldsymbol{\tau}|}{\rho_a}}.
$$

### High-Frequency Stress Tail

For accurate form stress, the manual recommends a high-frequency limit
of at least 2 Hz. The implementation appends a wind-speed-dependent
power-law stress tail from the highest resolved wavenumber to a
capillary wavenumber:

$$
p(U) = aU^2+bU+c,
$$

with

$$
a=0.000112,\qquad b=-0.01451,\qquad c=-1.0186.
$$

The tail factor is

$$
T
=
\frac{k_{cap}^{p+1}-k_{max}^{p+1}}
{k_{max}^{p}(p+1)},
\qquad
k_{cap}=10^3\,\mathrm{m^{-1}}.
$$

This tail is applied in the wind direction for atmospheric stress and in
the dissipative direction for ocean-top momentum flux.

## Ocean Momentum Fluxes

Momentum fluxes into the ocean are defined positive downward. Ocean-top
flux receives dissipation by breaking, turbulence, and viscosity, plus
tail and skin stress:

$$
\tau_x^{OT}
=
\rho_w g
\int_{-\pi}^{\pi}
\int_{k_{\min}}^{k_{\max}}
\frac{-S_{ds}-S_{dt}-S_{dv}}{c}
\cos\phi\,k\,dk\,d\phi
+ \tau_x(tail)
+ \tau_x(skin),
$$

$$
\tau_y^{OT}
=
\rho_w g
\int_{-\pi}^{\pi}
\int_{k_{\min}}^{k_{\max}}
\frac{-S_{ds}-S_{dt}-S_{dv}}{c}
\sin\phi\,k\,dk\,d\phi
+ \tau_y(tail)
+ \tau_y(skin).
$$

Ocean-bottom flux receives bottom friction and bottom percolation:

$$
\tau_x^{OB}
=
\rho_w g
\int_{-\pi}^{\pi}
\int_{k_{\min}}^{k_{\max}}
\frac{-S_{bf}-S_{bp}}{c}
\cos\phi\,k\,dk\,d\phi,
$$

$$
\tau_y^{OB}
=
\rho_w g
\int_{-\pi}^{\pi}
\int_{k_{\min}}^{k_{\max}}
\frac{-S_{bf}-S_{bp}}{c}
\sin\phi\,k\,dk\,d\phi.
$$

The source additionally diagnoses `taux_snl` and `tauy_snl`, the
momentum loss associated with nonlinear downshifting. `Snl` conserves
energy in the implemented transfer but does not conserve momentum
because phase speed changes between source and receiver bins.

## Stokes Drift

When `OUTPUT: stokes = .TRUE.`, the model reads user-provided depth
levels from the `STOKES` namelist and computes wave-induced Stokes drift
at those levels.

Depths in the namelist are positive downward. Internally they are stored
as $z_l \le 0$, with $z=0$ at the surface and $z=-d$ at the bottom.

For each requested level, UMWM computes

$$
\mathbf{u}_S(z_l)
=
\sum_{o,p}
K_{o,l}
E_{o,p}
\begin{bmatrix}
\cos\phi_p\\
\sin\phi_p
\end{bmatrix},
$$

where, for finite depth,

$$
K_{o,l}
=
\omega_o k_o^2
\frac{\cosh\left[2k_o(z_l+d)\right]}{\sinh^2(k_od)}
dk_o\,d\phi.
$$

For large arguments that would overflow hyperbolic functions, the code
uses the deep-water approximation

$$
K_{o,l}
=
2\omega_o k_o^2
\exp(2k_o z_l)
dk_o\,d\phi.
$$

If a requested level is below the local bottom, its Stokes contribution
is set to zero.

The output `d_stokes` is the e-folding depth of the Stokes drift
magnitude, found by linear interpolation to the depth where

$$
|\mathbf{u}_S(z)| = \frac{|\mathbf{u}_S(0)|}{e}.
$$

## Integrated Diagnostics

The following diagnostics are computed from the full spectrum unless
noted otherwise.

Significant wave height:

$$
H_s
=
4\sqrt{
\int\!\int E(k,\phi)k\,dk\,d\phi
}.
$$

Discrete implementation:

$$
H_s
=
4\sqrt{
\sum_{o,p}E_{o,p} k_o\,dk_o\,d\phi
}.
$$

Mean wave period:

$$
\overline{T}
=
\sqrt{
\frac{
\int\!\int E(k,\phi)k\,dk\,d\phi
}{
\int\!\int f^2E(k,\phi)k\,dk\,d\phi
}
}.
$$

Mean wavelength:

$$
\overline{L}
=
2\pi
\sqrt{
\frac{
\int\!\int E(k,\phi)k\,dk\,d\phi
}{
\int\!\int k^2E(k,\phi)k\,dk\,d\phi
}
}.
$$

Mean wave direction:

$$
\overline{\phi}
=
\operatorname{atan2}
\left(
\sum_p M_p\sin\phi_p,
\sum_p M_p\cos\phi_p
\right),
$$

where

$$
M_p = \sum_o E_{o,p}k_o\,dk_o.
$$

Mean-square slope:

$$
mss =
\sum_{o,p} E_{o,p}k_o^3\,dk_o\,d\phi.
$$

The dominant period, wavelength, and direction are taken from the
discrete spectral peak of

$$
E_{o,p}k_o\,dk_o.
$$

The wave momentum diagnostics are integrated over the prognostic range:

$$
M_x
=
\rho_w g
\sum_{o\le oc,p}
E_{o,p}
\frac{k_o\,dk_o}{c_o}
\cos\phi_p\,d\phi,
$$

$$
M_y
=
\rho_w g
\sum_{o\le oc,p}
E_{o,p}
\frac{k_o\,dk_o}{c_o}
\sin\phi_p\,d\phi.
$$

The radiation-like momentum flux components are

$$
C_{g}M_{xx}
=
\rho_w g
\sum_{o\le oc,p}
c_{g,o}E_{o,p}
\frac{k_o\,dk_o}{c_o}
\cos^2\phi_p\,d\phi,
$$

$$
C_{g}M_{xy}
=
\rho_w g
\sum_{o\le oc,p}
c_{g,o}E_{o,p}
\frac{k_o\,dk_o}{c_o}
\cos\phi_p\sin\phi_p\,d\phi,
$$

$$
C_{g}M_{yy}
=
\rho_w g
\sum_{o\le oc,p}
c_{g,o}E_{o,p}
\frac{k_o\,dk_o}{c_o}
\sin^2\phi_p\,d\phi.
$$

## Physics Namelist Controls

The main physics controls live in the `PHYSICS` namelist. Defaults below
are from `namelists/main.nml`.

| Parameter | Default | Units | Role |
| --- | ---: | --- | --- |
| `g` | `9.80665` | `m/s^2` | Gravitational acceleration. |
| `nu_air` | `1.56E-5` | `m^2/s` | Air kinematic viscosity for skin roughness. |
| `nu_water` | `0.90E-6` | `m^2/s` | Water kinematic viscosity for `Sdv`. |
| `sfct` | `0.07` | `N/m` | Surface tension in dispersion relation. |
| `kappa` | `0.4` | nondimensional | Von Karman constant. |
| `z` | `10.` | `m` | Input wind height. |
| `gustiness` | `0.0` | nondimensional | Random wind-component perturbation fraction; should be nonnegative and not larger than about `0.2`. |
| `dmin` | `10.` | `m` | Minimum depth applied to ocean cells. |
| `explim` | `0.9` | nondimensional | Maximum source exponent used to choose `dts`; `0.69` permits about 100 percent growth. |
| `sin_fac` | `0.11` | nondimensional | Wind-sea growth factor. |
| `sin_diss1` | `0.10` | nondimensional | Opposing-wind damping factor. |
| `sin_diss2` | `0.001` | nondimensional | Swell-overrunning-wind damping factor. |
| `sds_fac` | `42.` | nondimensional | Breaking dissipation factor $A_2$. |
| `sds_power` | `2.4` | nondimensional | Saturation-spectrum exponent $n$. |
| `mss_fac` | `360` | nondimensional | Mean-square-slope enhancement factor $A_3$. |
| `snl_fac` | `5.0` | nondimensional | Nonlinear downshifting factor $A_5$. |
| `sdt_fac` | `0.002` | nondimensional | Turbulent dissipation factor $A_4$. |
| `sbf_fac` | `0.003` | `m/s` | Bottom friction coefficient $G_f$. |
| `sbp_fac` | `0.003` | `m/s` | Bottom percolation coefficient $G_p$. |
| `fcut_fac` | `0.53` | nondimensional | Factor of the prognostic cutoff frequency $f_c$ = `fcut_fac`$\,g/U_{10}$, about 4 times the Pierson-Moskowitz peak frequency by default. |

Physics-related `DOMAIN` controls:

| Parameter | Role |
| --- | --- |
| `om` | Number of frequency/wavenumber bins. |
| `pm` | Number of direction bins. The manual requires divisibility by `4`; the CFL optimization is best when divisible by `8`. |
| `fmin`, `fmax` | Lowest and highest frequency bins. |
| `fprog` | Upper cap on the prognostic frequency range. The local cutoff is `min(fcut_fac*g/U10, fprog)`. |
| `dtg` | Global forcing/output step and maximum source-loop interval. |

Sea-ice controls:

| Parameter | Default | Role |
| --- | ---: | --- |
| `fice0` | `0.0` | Constant sea-ice fraction when sea ice is not read from input files. |
| `fice_lth` | `0.30` | Lower sea-ice fraction threshold for attenuation. |
| `fice_uth` | `0.75` | Upper sea-ice fraction threshold above which propagation/refraction zero the updated spectrum. |

`OUTPUT` namelist controls:

| Parameter | Default | Units | Role |
| --- | ---: | --- | --- |
| `outgrid` | `1` | hours | Gridded output interval. A value of `0` disables gridded output. |
| `outspec` | `0` | hours | Spectrum point output interval. A value of `0` disables spectrum output. |
| `outrst` | `6` | hours | Restart output interval. A value of `0` disables restart output. |
| `xpl` | `20` | grid index | Grid-cell x index used for stdout diagnostic reporting. |
| `ypl` | `6` | grid index | Grid-cell y index used for stdout diagnostic reporting. |
| `stokes` | `.t.` | logical | Enables Stokes drift calculation and output; when true, the `STOKES` namelist is read. |

`STOKES` namelist controls:

| Parameter | Default | Units | Role |
| --- | --- | --- | --- |
| `depths` | `0.1 0.2 ... 100` | m | Positive-down depth levels where `u_stokes`, `v_stokes`, and `d_stokes` diagnostics are computed. |

## Output Fields Relevant to Physics

Gridded output includes forcing, integrated wave quantities, stress, and
diagnostics:

| Field | Meaning |
| --- | --- |
| `wspd`, `wdir` | Wind speed and mathematical wind direction. |
| `uc`, `vc` | Ocean current components. |
| `rhoa`, `rhow` | Air and water density. |
| `fice` | Sea-ice fraction. |
| `taux_form`, `tauy_form` | Wind form drag components. |
| `taux_skin`, `tauy_skin` | Skin drag components. |
| `taux_ocn`, `tauy_ocn` | Momentum flux into ocean top. |
| `taux_bot`, `tauy_bot` | Momentum flux into ocean bottom. |
| `taux_snl`, `tauy_snl` | Momentum loss due to nonlinear downshifting. |
| `cd` | Total atmospheric drag coefficient. |
| `ust` | Air-side friction velocity. |
| `shelt` | Variable sheltering coefficient. |
| `swh`, `mwp`, `mwl`, `mwd` | Significant height, mean period, mean wavelength, mean direction. |
| `mss` | Mean-square slope. |
| `dwp`, `dwl`, `dwd` | Dominant period, wavelength, and direction from the discrete spectral peak. |
| `momx`, `momy` | Wave momentum components. |
| `cgmxx`, `cgmxy`, `cgmyy` | Group-velocity-weighted momentum flux components. |
| `epsx_atm`, `epsy_atm` | Wave energy growth flux from air. |
| `epsx_ocn`, `epsy_ocn` | Wave energy dissipation flux to ocean. |
| `physics_time_step` | Local source-limited physics timestep diagnostic. |
| `u_stokes`, `v_stokes`, `d_stokes` | Stokes drift components and e-folding depth, when enabled. |

Spectrum point output includes:

| Field | Meaning |
| --- | --- |
| `F` | Surface elevation variance spectrum. |
| `Sin` | Wind input/export implementation rate. |
| `Sds` | Wave breaking implementation decay rate. |
| `Sdt` | Turbulent implementation decay rate. |
| `Sdv` | Viscous implementation decay rate. |
| `Sbf` | Combined bottom friction and bottom percolation implementation decay rate. |
| `Snl` | Nonlinear wave-wave interaction absolute source/sink tendency. |

The current code writes the implementation arrays directly in spectrum
point output. Thus `Sin`, `Sds`, `Sdt`, `Sdv`, and `Sbf` are rates that
must be multiplied by `F` to obtain variance-spectrum tendencies, while
`Snl` is already an absolute tendency.

## Tests

The focused unit tests live in `tests/` and are built by
`tests/Makefile`. They can be run either from the top level

```sh
make test
```

or directly from the test directory

```sh
make --directory=tests test
```

The current focused tests are:

| Test | Coverage |
| --- | --- |
| `test_dispersion` | Wavenumber/frequency round-trip accuracy and positive finite group speed for the capillary-gravity dispersion module. |
| `test_spectrum` | Spectrum constructor metadata, endpoint frequencies, legacy log-space frequency equivalence, direction-grid equivalence, `dlnf`, and `dth`. |
| `test_source_functions` | Mean-squared slope moments against the direct sum over directions, the vectorized power function against the intrinsic, and the wind input against a step-by-step evaluation. |
| `test_domain` | One-hour integration of the idealized hurricane case by `domain_type` with gridded and spectrum output and Stokes drift, compared with the mean and maximum of output fields, ocean stress, and Stokes drift from UMWM before the domain refactor. |

The Python package has its own test suite in `python/tests`, run with
`python -m pytest python/tests`, which includes the same idealized
hurricane case driven from Python and compared with the same reference
values, and `python/tests/mpi_smoke.py` for a library built with MPI.

`make clean` at the top level recurses into `tests/` and removes
test-local binaries, object files, and module files.

## References

The physics basis and manual equations are from UMWM manual version 2,
sections 1 through 4 and 9. The cited references are:

- Donelan, M. A., and W. J. Plant, 2009: A threshold for wind-wave growth. *Journal of Geophysical Research*, 114, C07012. https://doi.org/10.1029/2008JC005238.
- Donelan, M. A., M. Curcic, S. S. Chen, and A. K. Magnusson, 2012: Modeling waves and wind stress. *Journal of Geophysical Research*, 117, C00J23. https://doi.org/10.1029/2011JC007787.
- Jeffreys, H., 1924: On the formation of waves by wind. *Proceedings of the Royal Society A*, 107, 189-206.
- Jeffreys, H., 1925: On the formation of waves by wind, II. *Proceedings of the Royal Society A*, 110, 341-347.
- Komen, G. J., L. Cavaleri, M. Donelan, K. Hasselmann, S. Hasselmann, and P. A. E. M. Janssen, 1994: *Dynamics and Modelling of Ocean Waves*. Cambridge University Press, 532 pp.
- Pierson, W. L., and L. Moskowitz, 1964: A proposed spectral form for fully developed wind seas based on the similarity theory of S. A. Kitaigorodskii. *Journal of Geophysical Research*, 69, 5181-5190.
- Shemdin, O. H., K. Hasselmann, S. V. Hsiao, and K. Herterich, 1978: Non-linear and linear bottom interaction effects in shallow water. In A. Favre and K. Hasselmann (eds.), *Turbulent Fluxes Through the Sea Surface, Wave Dynamics, and Prediction*, Plenum, New York, 347-372.
