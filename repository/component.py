"""Repository file to component with create, delete and get functions
"""
import asyncio
from database.tortoise.models import Component
from database.excell.component_data import COMPONENTS


async def create_components():
    """A function to create component in datbase
    """
    components_data = list(COMPONENTS.values())
    components = await asyncio.gather(*[Component.create(**component) for component in components_data])
    for component in components:
        await component.save()


async def delete_all_components():
    """A function to delete all fluids in datbase
    """
    await Component.filter().delete()


async def get_all_components():
    """A function to get all fluids in datbase
    """
    return await Component.all()