"""file that contains seed function"""
import asyncio
from database.tortoise.models import Component, Fluid, PengRobinsonParam, ComponentFluid
from repository.fluid import get_all_fluids, create_fluids, delete_all_fluids
from repository.component import create_components, delete_all_components, get_all_components
from repository.peng_robinson_params import create_peng_robinson_params, delete_all_peng_robinson_params, get_all_peng_robinson_params
from repository.composition import create_compositions, delete_all_compositions, get_all_compositions


async def seed():
    """function to generate data for project
    """
    await asyncio.gather(delete_all_fluids(), delete_all_components(), delete_all_peng_robinson_params(), delete_all_compositions())
    await asyncio.gather(create_components(), create_fluids())
    await asyncio.gather(create_peng_robinson_params(), create_compositions())

    fluids = await get_all_fluids()
    for fluid in fluids:
        print(fluid.name)
        for field in Fluid._meta.fields_map:
            print(f"  {field}: {getattr(fluid, field)}")

    components = await get_all_components()
    for component in components:
        print(component.name)
        for field in Component._meta.fields_map:
            print(f"  {field}: {getattr(component, field)}")

    peng_robinson_params = await get_all_peng_robinson_params()
    for peng_robinson_param in peng_robinson_params:
        for field in PengRobinsonParam._meta.fields_map:
            print(f"  {field}: {getattr(peng_robinson_param, field)}")

    compositions = await get_all_compositions()
    for composition in compositions:
        for field in ComponentFluid._meta.fields_map:
            print(f"  {field}: {getattr(composition, field)}")
