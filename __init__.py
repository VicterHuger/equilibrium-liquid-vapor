import asyncio
from dotenv import load_dotenv
from tortoise import Tortoise, run_async
from database.tortoise.tortoise_config import TORTOISE
from database.tortoise.models import Fluid
from database.excell.data import create_fluid_data

load_dotenv()


async def create_tables():
    await Tortoise.init(
        config=TORTOISE,
        modules={"models": ["database.tortoise.models"]}
    )
    await Tortoise.generate_schemas()


async def create_fluids():
    fluids_data = [create_fluid_data('Fluid1')]
    fluids_data.append(create_fluid_data('Fluid2'))
    fluids = await asyncio.gather(*[Fluid.create(**fluid) for fluid in fluids_data])
    for fluid in fluids:
        await fluid.save()
        print(f"Fluid {fluid.id}:")
        for field in Fluid._meta.fields_map:
            print(f"  {field}: {getattr(fluid, field)}")


async def delete_all_fluids():
    await Fluid.filter().delete()


async def init():
    await create_tables()
    print('Tables Created')
    await delete_all_fluids()
    await create_fluids()
    fluids = await Fluid.all()
    for fluid in fluids:
        print(f"Fluid {fluid.id}:")
        for field in Fluid._meta.fields_map:
            print(f"  {field}: {getattr(fluid, field)}")

    await Tortoise.close_connections()

run_async(init())
