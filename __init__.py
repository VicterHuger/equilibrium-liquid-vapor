"""file to inialize project"""
from dotenv import load_dotenv
from tortoise import Tortoise, run_async
from database.tortoise.tortoise_config import TORTOISE
from database.tortoise.models import Fluid
from repository.fluid import get_all_fluids, create_fluids, delete_all_fluids

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
    await delete_all_fluids()
    await create_fluids()
    fluids = await get_all_fluids()
    for fluid in fluids:
        print(fluid.name)
        for field in Fluid._meta.fields_map:
            print(f"  {field}: {getattr(fluid, field)}")

    await Tortoise.close_connections()

run_async(init())
