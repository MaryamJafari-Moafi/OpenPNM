import OpenPNM.openpnm.models.physics as mods

Dless_numbers = {
    'pore.reynolds_number': {
        'model': mods.hydraulics.reynolds_number,
    },
    'pore.peclet_number': {
        'model': mods.hydraulics.peclet_number,
    },
    'pore.prandtl_number': {
        'model': mods.hydraulics.prandtl_number,
    },
    'pore.nusselt_number': {
        'model': mods.hydraulics.nusselt_number,
    },
}