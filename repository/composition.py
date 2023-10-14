"""Repository file for composition creating, deleting and getting operations
"""
from database.tortoise.models import ComponentFluid
from database.seed_data.composition_data import COMPOSITIONS_VALUES
from repository.fluid import get_fluid
from repository.component import get_component


async def create_compositions():
    """A function to create composition in database
    """

    compositions_data = COMPOSITIONS_VALUES.items()
    for fluid_name, compositions_dict in compositions_data:
        fluid_obj = await get_fluid(fluid_name)
        for component_name, composition_value in compositions_dict.items():
            component_obj = await get_component(component_name)
            composition = ComponentFluid(
                component=component_obj, fluid=fluid_obj, composition=composition_value)
            await composition.save()


async def delete_all_compositions():
    """A function to delete all compositions in database
    """
    await ComponentFluid.filter().delete()


async def get_all_compositions():
    """A function to get all compositions in database
    """
    return await ComponentFluid.all()
