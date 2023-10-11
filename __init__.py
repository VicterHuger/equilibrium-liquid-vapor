"""file to inialize project"""
import asyncio
from dotenv import load_dotenv
from tortoise import Tortoise, run_async
from database.tortoise.tortoise_config import TORTOISE
from database.tortoise.models import Component, Fluid
from repository.fluid import get_all_fluids, create_fluids, delete_all_fluids
from repository.component import create_components, delete_all_components, get_all_components

load_dotenv()


async def create_tables():
    """function to create tortoise tables in database
    """
    await Tortoise.init(
        config=TORTOISE,
        modules={"models": ["database.tortoise.models"]}
    )
    await Tortoise.generate_schemas()


async def init():
    """function to inializr project
    """
    await create_tables()
    print('Tables Created')
    await asyncio.gather(delete_all_fluids(), delete_all_components())
    await asyncio.gather(create_components(), create_fluids())

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

    await Tortoise.close_connections()

run_async(init())
