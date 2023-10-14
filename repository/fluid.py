"""Repository file to fluid with create, delete and get functions
"""
import asyncio
from database.tortoise.models import Fluid
from database.seed_data.fluid_data import create_fluid_data


async def create_fluids():
    """A function to create fluid in datbase
    """
    fluids_data = [create_fluid_data('Fluid1')]
    fluids_data.append(create_fluid_data('Fluid2'))
    fluids_data.append(create_fluid_data('Fluid3'))
    fluids_data.append(create_fluid_data('Fluid4'))
    fluids_data.append(create_fluid_data('Fluid5'))
    fluids = await asyncio.gather(*[Fluid.create(**fluid) for fluid in fluids_data])
    for fluid in fluids:
        await fluid.save()


async def delete_all_fluids():
    """A function to delete all fluids in datbase
    """
    await Fluid.filter().delete()


async def get_all_fluids():
    """A function to get all fluids in datbase
    """
    return await Fluid.all()
