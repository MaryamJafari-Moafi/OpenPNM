import numpy as _np
from openpnm.models import _doctxt

__all__ = [
    "reynolds_number",
    "peclet_number",
    "prandtl_number",
    "nusselt_number",
]


@_doctxt
def reynolds_number(
    phase,
    pore_velocity="pore.velocity",
    pore_density="pore.density",
    pore_viscosity="pore.viscosity",
    pore_volume="pore.volume",
    pore_surface_area="pore.surface_area",
    length="throat.diameter",
):
    r"""
    Calculates the Reynolds number for fluid in each pore in a porous medium.

    The Reynolds number is calculated as:

    .. math::

        Re = \frac{\rho \cdot v \cdot d_h}{\mu}

    where:
    - :math:`\rho` is the fluid density,
    - :math:`v` is the velocity of the fluid,
    - :math:`d_h` is the hydraulic diameter (computed as 4 * pore.volume / pore.surface_area),
    - :math:`\mu` is the dynamic viscosity.

    Parameters
    ----------
    %(phase)s
    pore_velocity : str
        %(dict_blurb)s velocity of the fluid in the pores.
    pore_density : str
        %(dict_blurb)s density of the fluid in the pores.
    pore_viscosity : str
        %(dict_blurb)s dynamic viscosity of the fluid in the pores.
    pore_volume : str
        %(dict_blurb)s volume of the pores.
    pore_surface_area : str
        %(dict_blurb)s surface area of the pores.

    Returns
    -------
    %(return_arr)s Reynolds number for each pore.
    """

    v = phase[pore_velocity]
    rho = phase[pore_density]
    mu = phase[pore_viscosity]
    sa = phase[pore_surface_area]
    vol = phase[pore_volume]
    h_d = 4 * vol / sa   #hydraulic_diameter
    return rho * v * h_d / mu


@_doctxt
def peclet_number(
    phase,
    pore_velocity="pore.velocity",
    pore_density="pore.density",
    pore_volume="pore.volume",
    pore_surface_area="pore.surface_area",
    pore_specific_heat="pore.heat_capacity",
    thermal_conductivity="pore.thermal_conductivity",
):
    r"""
    Calculates the Peclet number for heat transfer in a porous medium.

    The Peclet number is calculated as:

    .. math::

        Pe = \frac{\rho \cdot v \cdot d_h \cdot C_p}{k}

    where:
    - :math:`\rho` is the fluid density,
    - :math:`v` is the velocity of the fluid,
    - :math:`d_h` is the hydraulic diameter (computed as 4 * volume / surface area),
    - :math:`C_p` is the specific heat capacity of the fluid,
    - :math:`k` is the thermal conductivity of the fluid.

    Parameters
    ----------
    %(phase)s
    pore_velocity : str
        %(dict_blurb)s velocity of the fluid in the pores.
    pore_density : str
        %(dict_blurb)s density of the fluid in the pores.
    pore_volume : str
        %(dict_blurb)s volume of the pores.
    pore_surface_area : str
        %(dict_blurb)s surface area of the pores.
    pore_specific_heat : str
        %(dict_blurb)s specific heat capacity of the fluid in the pores.
    thermal_conductivity : str
        %(dict_blurb)s thermal conductivity of the fluid in the pores.

    Returns
    -------
    %(return_arr)s Peclet number for each pore.
    """
        
    v = phase[pore_velocity]
    rho = phase[pore_density]
    sa = phase[pore_surface_area]
    vol = phase[pore_volume]
    cp = phase[pore_specific_heat]
    k = phase[thermal_conductivity]
    h_d = 4 * vol / sa   #hydraulic_diameter
    return rho * v * h_d * cp / k


@_doctxt
def prandtl_number(
    phase,
    pore_viscosity="pore.viscosity",
    pore_specific_heat="pore.heat_capacity",
    thermal_conductivity="pore.thermal_conductivity"
):
    r"""
    Calculates the Prandtl number for a fluid.

    The Prandtl number is calculated as:

    .. math::

        Pr = \frac{\mu \cdot C_p}{k}

    where:
    - :math:`\mu` is the dynamic viscosity,
    - :math:`C_p` is the specific heat capacity of the fluid,
    - :math:`k` is the thermal conductivity.

    Parameters
    ----------
    %(phase)s
    pore_viscosity : str
        %(dict_blurb)s dynamic viscosity of the fluid in the pores.
    pore_specific_heat : str
        %(dict_blurb)s specific heat capacity of the fluid in the pores.
    thermal_conductivity : str
        %(dict_blurb)s thermal conductivity of the fluid in the pores.

    Returns
    -------
    %(return_arr)s Prandtl number for each pore.
    """
    mu = phase[pore_viscosity]
    cp = phase[pore_specific_heat]
    k = phase[thermal_conductivity]
    return mu * cp / k


@_doctxt
def nusselt_number(
    phase,
    reynolds="pore.reynolds_number",
    prandtl="pore.prandtl_number"
):
    r"""
    Calculates the Nusselt number using an empirical correlation for water, air and oil [1].

    The correlation is valid for:
    - Reynolds numbers in the range :math:`3.5 \leq Re \leq 76,000`
    - Prandtl numbers in the range :math:`0.71 \leq Pr \leq 380`

    The Nusselt number is calculated as:

    .. math::

        Nu = 2 + (0.4 \cdot Re^{0.5} + 0.06 \cdot Re^{\frac{2}{3}}) \cdot Pr^{0.4}

    Parameters
    ----------
    %(phase)s
    reynolds : str
        %(dict_blurb)s Reynolds number for each pore.
    prandtl : str
        %(dict_blurb)s Prandtl number for each pore.

    Returns
    -------
    %(return_arr)s Nusselt number for each pore.

    References
    ----------
    [1] Whitaker S. Forced convection heat transfer correlations for flow in
    pipes, past flat plates, single. AIChE J. 1972;18(2):361-371.
    """
    Re = phase[reynolds]
    Pr = phase[prandtl]
    return 2 + (0.4 * Re**0.5 + 0.06 * Re**(2/3)) * Pr**0.4  # correlation
